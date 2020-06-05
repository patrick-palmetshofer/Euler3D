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

void Boundary::setFront()
{
	dim = 2;
	dir = 1;
}

void Boundary::setBack()
{
	dim = 2;
	dir = -1;
}

void Boundary::setDimDir(int dim, int dir)
{
	this->dim = dim;
	this->dir = dir;
}

void Boundary::init()
{
	if (cell_inds.empty())
	{
		int dirind = 0;
		
		//For dim = 0; offdims are 1,2 and so on
		int nondim1 = 0, nondim2 = 1;
		switch (dim)
		{
		case 0:
			nondim1++;
		case 1:
			nondim2++;
		}

		//Get directional index
		if (dir == -1)
			dirind = grid->getnComponentCells(dim);

		cell_inds.reserve(grid->getnComponentCells(nondim1)*grid->getnComponentCells(nondim2));
		for (int i1 = 0; i1 < grid->getnComponentCells(nondim1); i1++)
		{
			for (int i2 = 0; i2 < grid->getnComponentCells(nondim2); i2++)
			{
				IndArray temp(0,0,0);
				temp[dim] = dirind;
				temp[nondim1] = i1;
				temp[nondim2] = i2;
				cell_inds.push_back(temp);
			}
		}
	}
	dir_arr.fill(0);
	dir_arr[dim] = dir;
}

void Boundary::apply()
{
	std::pair<StateVector, StateVector> leftrightstates;

	std::vector<StateVector> fl_test;
	fl_test.reserve((cell_inds.size()));

	for (auto &ind : cell_inds)
	{
		if (dir == 1)
		{
			leftrightstates.second = (*conservative)(ind[0])(ind[1])(ind[2]);
			leftrightstates.first = getGhostState(ind);
		}
		else if (dir == -1)
		{
			IndArray locind = ind + dir_arr;
			leftrightstates.first = (*conservative)(locind[0])(locind[1])(locind[2]);
			leftrightstates.second = getGhostState(locind);
		}

		leftrightstates.first = Euler::swap(leftrightstates.first, dim + 1);
		leftrightstates.second = Euler::swap(leftrightstates.second, dim + 1);

		DirVector n0 = grid->getnVec(ind, dim);
		DirVector n = Euler::swap(n0, 0, dim);

		StateVector flswapped = flux->calcFlux(leftrightstates, n);
		StateVector fl = Euler::swap(flswapped, dim + 1);
		if (Euler::checkNaN(&fl))
			throw;
		flux->setBoundaryFlux(ind, fl);
	}
	auto a = fl_test.size();
}

StateVector SupersonicInlet::getGhostState(IndArray& ind)
{
	return in_state;
}

StateVector SlipWall::getGhostState(IndArray& ind)
{
	DirVector n = grid->getnVec(ind, dim);

	double unorm = (*conservative)(ind[0])(ind[2])(ind[2]).matrix().segment(1,3).dot(n);

	ghost_state = (*conservative)(ind[0])(ind[2])(ind[2]);
	ghost_state.segment(1, 3) -= 2 * unorm*n.array();

	return ghost_state;
}

StateVector SupersonicOutlet::getGhostState(IndArray &ind)
{
	DirVector n = grid->getnVec(ind, dim);
	ghost_state = (*conservative)(ind[0])(ind[2])(ind[2]);

	return ghost_state;
}

StateVector BoundaryLODI::getGhostState(IndArray &ind)
{
	DirVector n = grid->getnVec(ind, dim);

	IndArray indpos = ind + dir_arr;

	StateVector &c = (*conservative)(indpos[0])(indpos[1])(indpos[2]);
	StateVector &cneg = (*conservative)(ind[0])(ind[1])(ind[2]);

	/*	cpos = c + (c - cneg);*/
	double u = c[2] / c[0];
	double sound = fluid->calcSoundSpeedcons(c);
	double rho = c[0];

	double dx = grid->getPoint(indpos)[dim]- grid->getPoint(ind)[dim];
	if (dx == 0)
		dx = 1;

	double drho = c[0] - cneg[0];
	double dp = fluid->calcPcons(c) - fluid->calcPcons(cneg);
	double du = (c[1] / c[0] - cneg[1] / cneg[0]);

	calcWaveStrengths(u,sound,rho,dx,drho,dp,du,c,cneg);

	double p = fluid->calcPcons(c);

	double ppos = p + dx * 0.5*(L5 / (u + sound) + L1 / (u - sound));
	StateVector prim = fluid->cons2prim(c);
	StateVector primpos;
	primpos[0] = prim[0] + dx * (L2 / u + 0.5*(L5 / (u + sound) + L1 / (u - sound))) / (sound*sound);
	primpos[1] = prim[1] + dx * (L5 / (u + sound) - L1 / (u - sound)) / (2 * sound*rho);
	primpos[2] = prim[2] + dx * L3 / u;
	primpos[3] = prim[3] + dx * L4 / u;
	primpos[4] = ppos / ((fluid->getGamma() - 1)*primpos[0]) + 0.5*primpos.segment(1, 3).matrix().squaredNorm();

	ghost_state = fluid->prim2cons(primpos);

	if (u == 0 || u + sound == 0 || u - sound == 0) {
		ghost_state = c + (c - cneg);
	}

	if (Euler::checkNaN(&ghost_state))
		throw;
	return ghost_state;
}

void Outlet::calcWaveStrengths(double u, double sound, double rho, double dx, double drho, double dp, double du, StateVector& c, StateVector& cneg)
{
	L1 = (u - sound)*(dp - rho * sound*du) / dx;
	L2 = u * (sound*sound*drho - dp) / dx;
	L3 = u * (c[2] / c[0] - cneg[2] / cneg[0]) / dx;
	L4 = u * (c[3] / c[0] - cneg[3] / cneg[0]) / dx;
	L5 = (u + sound)*(dp + rho * sound*du) / dx;

	double p = fluid->calcPcons(c);
	double Ma = u / sound;
	if (Ma*Ma < 1)
	{
		double K = 0.58*(1 - Ma * Ma)*sound / grid->getMaxDistance(dim);
		L1 = K * (p - p_infty);
	}
}
