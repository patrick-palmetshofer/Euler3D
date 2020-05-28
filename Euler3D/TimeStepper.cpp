#include "TimeStepper.h"



TimeStepper::TimeStepper()
{
}


TimeStepper::~TimeStepper()
{
}

//Calculates maximum time step for a single cell (does not apply maxCFL)
double TimeStepper::calcTimeStep(int i, int j, StateMatrix2D * conservative)
{
	p = fluid->cons2prim((*conservative)(i,j));
	double c = fluid->calcSoundSpeedcons((*conservative)(i,j));// std::sqrt(gamma*(gamma - 1)*(p[3] - 0.5*(p[1] * p[1] + p[2] * p[2])));

	double Uxi_velnorm = grid->getnXiXs(i,j) * p[1] + grid->getnXiYs(i, j) * p[2];
	double Ueta_velnorm = grid->getnEtaXs(i, j) * p[1] + grid->getnEtaYs(i, j) * p[2];

	double r_xi = std::abs(Uxi_velnorm) + c;
	double r_eta = std::abs(Ueta_velnorm) + c;

	//	double dtlocal = std::min((points(i+1,j)[0] - points(i,j)[0]) / r_xi, (points(i,j+1)[1] - points(i,j)[1]) / r_eta);
	double dtlocal = std::min(grid->getSxi(i,j) / r_xi, grid->getSeta(i, j) / r_eta);

	return dtlocal;
}

//Calculates maximum time step for a all cells
double TimeStepper::calcTimeStep(StateMatrix2D * conservative)
{
	double dt = 1e19;
	int xi_size = conservative->rows();
	int eta_size = (*conservative).cols();
	for (int i = 0; i < xi_size; i++)
	{
		for (int j = 0; j < eta_size; j++)
		{
			//double dtlocal = std::min(1 / r_xi, 1 / r_eta);
			double dtlocal = calcTimeStep(i, j,conservative);
			if (dtlocal < dt)
				dt = dtlocal;
		}
	}
	return dt;
}