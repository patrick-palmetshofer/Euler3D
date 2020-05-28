#include "Boundary.h"



Boundary::Boundary()
{
}


Boundary::~Boundary()
{
}

void Boundary::setBottom()
{
	dim = 1;
	dir = 1;
}

void Boundary::setTop()
{
	dim = 1;
	dir = -1;
}

void Boundary::setRight()
{
	dim = 0;
	dir = -1;
}

void Boundary::setLeft()
{
	dim = 0;
	dir = 1;
}

void Boundary::init()
{
	if (cell_inds.empty())
	{
		if (dim == 0)
		{
			int dirind = 0;
			if (dir == -1)
				dirind = grid->getnxiCells();
			for (int i = 0; i < grid->getnetaCells(); i++)
				cell_inds.push_back({ dirind, i });
		}
		else if (dim == 1)
		{
			int dirind = 0;
			if (dir == -1)
				dirind = grid->getnetaCells();
			for (int i = 0; i < grid->getnxiCells(); i++)
				cell_inds.push_back({ i, dirind });
		}
		else
			throw;
	}
}

void Boundary::apply()
{
	std::pair<StateVector2D, StateVector2D> leftrightstates;
	int physdim1 = 0, physdim2 = 1;
	if (dim == 1) 
	{
		physdim1++; physdim2--;
	}

	for (auto &inds : cell_inds)
	{
		if (dir == 1)
		{
			leftrightstates.second = (*conservative)(inds[0],inds[1]);
			leftrightstates.first = getGhostState(inds);
		}
		else if (dir == -1)
		{
			std::array<int,2> locinds = { inds[0] - physdim2, inds[1] - physdim1 };
			leftrightstates.first = (*conservative)(locinds[0],locinds[1]);
			leftrightstates.second = getGhostState(locinds);
		}

		leftrightstates.first = Euler::swap(leftrightstates.first, dim + 1);
		leftrightstates.second = Euler::swap(leftrightstates.second, dim + 1);

		StateVector2D fl = flux->calcFlux(leftrightstates, grid->getnComponent(inds[0], inds[1], dim, physdim1), grid->getnComponent(inds[0], inds[1], dim, physdim2));
		flux->setBoundaryFlux(inds[0], inds[1], Euler::swap(fl,dim+1));
	}
}

StateVector2D SupersonicInlet::getGhostState(std::array<int, 2> ind)
{
	return in_state;
}

StateVector2D SlipWall::getGhostState(std::array<int, 2> ind)
{
	int i = ind[0];
	int j = ind[1];
	int ineg = i, jneg = j;
	double nx = 0, ny = 0;
	switch (dim)
	{
	case 0:
		ineg += dir;
		nx = grid->getnXiXs(i, j);
		ny = grid->getnXiYs(i, j);
		break;
	case 1:
		jneg += dir;
		nx = grid->getnEtaXs(i, j);
		ny = grid->getnEtaYs(i, j);
		break;
	}

	double unorm = (*conservative)(i,j)[1] * nx + (*conservative)(i,j)[2]* ny;
	ghost_state[0] = (*conservative)(i,j)[0];
	ghost_state[1] = (*conservative)(i,j)[1] - 2 * unorm*nx;
	ghost_state[2] = (*conservative)(i,j)[2] - 2 * unorm*ny;
	ghost_state[3] = (*conservative)(i,j)[3];

	return ghost_state;
}

StateVector2D SupersonicOutlet::getGhostState(std::array<int, 2> ind)
{
	int i = ind[0];
	int j = ind[1];
	int ineg = i, jneg = j;
	double nx = 0, ny = 0;
	switch (dim)
	{
	case 0:
		ineg += dir;
		nx = grid->getnXiXs(i, j);
		ny = grid->getnXiYs(i, j);
		break;
	case 1:
		jneg += dir;
		nx = grid->getnEtaXs(i, j);
		ny = grid->getnEtaYs(i, j);
		break;
	}
	ghost_state = (*conservative)(i,j);
	
	return ghost_state;
}

StateVector2D BoundaryLODI::getGhostState(std::array<int, 2> ind)
{
	int i = ind[0];
	int j = ind[1];
	int ineg = i, jneg = j;
	double nx = 0, ny = 0;
	switch (dim)
	{
	case 0:
		ineg += dir;
		nx = grid->getnXiXs(i, j);
		ny = grid->getnXiYs(i, j);
		break;
	case 1:
		jneg += dir;
		nx = grid->getnEtaXs(i, j);
		ny = grid->getnEtaYs(i, j);
		break;
	}

	StateVector2D &c = (*conservative)(i,j);
	StateVector2D &cneg = (*conservative)(ineg,jneg);

	/*	cpos = c + (c - cneg);*/
	double u = c[2] / c[0];
	double sound = fluid->calcSoundSpeedcons(c);
	double rho = c[0];

	double dx = (grid->getPoint(i,j)[dim] - grid->getPoint(ineg, jneg)[dim]);

	double drho = c[0] - cneg[0];
	double dp = fluid->calcPcons(c) - fluid->calcPcons(cneg);
	double du = (c[1] / c[0] - cneg[1] / cneg[0]);

	calcWaveStrengths(u,sound,rho,dx,drho,dp,du,c,cneg);

	double p = fluid->calcPcons(c);

	double ppos = p + dx * 0.5*(L5 / (u + sound) + L1 / (u - sound));
	StateVector2D prim = fluid->cons2prim(c);
	StateVector2D primpos;
	primpos[0] = prim[0] + dx * (L2 / u + 0.5*(L5 / (u + sound) + L1 / (u - sound))) / (sound*sound);
	primpos[1] = prim[1] + dx * (L5 / (u + sound) - L1 / (u - sound)) / (2 * sound*rho);
	primpos[2] = prim[2] + dx * L3 / u;
	primpos[3] = ppos / ((fluid->getGamma() - 1)*primpos[0]) + 0.5*(primpos[1] * primpos[1] + primpos[2] * primpos[2]);

	ghost_state = fluid->prim2cons(primpos);

	if (u == 0 || u + sound == 0 || u - sound == 0) {
		ghost_state = c + (c - cneg);
	}

	Euler::checkNaN(&ghost_state);
	return ghost_state;
}

void Outlet::calcWaveStrengths(double u, double sound, double rho, double dx, double drho, double dp, double du, StateVector2D& c, StateVector2D& cneg)
{
	L1 = (u - sound)*(dp - rho * sound*du) / dx;
	L2 = u * (sound*sound*drho - dp) / dx;
	L3 = u * (c[1] / c[0] - cneg[1] / cneg[0]) / dx;
	L5 = (u + sound)*(dp + rho * sound*du) / dx;

	double p = fluid->calcPcons(c);
	double Ma = u / sound;
	if (Ma*Ma < 1)
	{
		double K = 0.58*(1 - Ma * Ma)*sound / grid->getMaxDistance(dim);
		L1 = K * (p - p_infty);
	}
}
