#include "TimeStepper.h"



TimeStepper::TimeStepper()
{
}


TimeStepper::~TimeStepper()
{
}

//Calculates maximum time step for a single cell (does not apply maxCFL)
double TimeStepper::calcTimeStep(int i, int j, int k, StateTensor * conservative)
{
	p = fluid->cons2prim((*conservative)(i)(j)(k));
	double c = fluid->calcSoundSpeedcons((*conservative)(i)(j)(k));// std::sqrt(gamma*(gamma - 1)*(p[3] - 0.5*(p[1] * p[1] + p[2] * p[2])));
	double dtlocal = 1e9;
	for (int dim = 0; dim < DirVector().size(); ++dim)
	{
		IndArray ind(i, j, k);
		double Ucomp_velnorm = p.segment(1,3).matrix().dot(grid->getnVec(ind, dim));

		double r = std::abs(Ucomp_velnorm) + c;

		//	double dtlocal = std::min((points(i+1,j)[0] - points(i)(j)(k)[0]) / r_xi, (points(i,j+1)[1] - points(i)(j)(k)[1]) / r_eta);
		dtlocal = std::min(grid->getFaceArea(i, j, k, dim) / r, dtlocal);
	}
	return dtlocal;
}

//Calculates maximum time step for a all cells
double TimeStepper::calcTimeStep(StateTensor * conservative)
{
	double dt = 1e19;
	int xi_size = conservative->rows();
	int eta_size = (*conservative).cols();
	for (int i = 0; i < conservative->size(); i++)
	{
		for (int j = 0; j < (*conservative)(0).size(); j++)
		{
			for (int k = 0; k < (*conservative)(0)(0).size(); k++)
			{
				//double dtlocal = std::min(1 / r_xi, 1 / r_eta);
				double dtlocal = calcTimeStep(i, j, k, conservative);
				if (dtlocal < dt)
					dt = dtlocal;
			}
		}
	}
	return dt;
}