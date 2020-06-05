#include "TimeStepperGlobal.h"

TimeStepperGlobal::TimeStepperGlobal()
{
}


TimeStepperGlobal::~TimeStepperGlobal()
{
}


void TimeStepperGlobal::execute(StateTensor * conservative, Eigen::Array<StateTensor *, 3, 1> fluxes)
{
	//Calculate time step with CFL
	dt = maxCFL * calcTimeStep(conservative);

	for (int i = 0; i < conservative->size(); i++)
	{
		for (int j = 0; j < (*conservative)(0).size(); j++)
		{
			for (int k = 0; k < (*conservative)(0)(0).size(); k++)
			{
				IndArray ind(i, j, k);
				diffs.fill(0);
				for (int dim = 0; dim < diffs.size(); ++dim)
				{
					IndArray indplus = ind;
					indplus[dim]++;
					diffs  += (*fluxes(dim))(indplus[0])(indplus[1])(indplus[2]) * grid->getFaceArea(indplus, dim) - (*fluxes(dim))(ind[0])(ind[1])(ind[2])* grid->getFaceArea(ind, dim);
				}
				(*conservative)(i)(j)(k) -= dt / grid->getVolume(i, j, k) * diffs;
			}
		}
	}

	time += dt;
}