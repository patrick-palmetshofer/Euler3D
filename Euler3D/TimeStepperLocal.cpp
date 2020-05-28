#include "TimeStepperLocal.h"



TimeStepperLocal::TimeStepperLocal()
{
}


TimeStepperLocal::~TimeStepperLocal()
{
}

void TimeStepperLocal::execute(StateMatrix2D * conservative, StateMatrix2D *xi_fluxes, StateMatrix2D *eta_fluxes)
{
	int xi_size = conservative->rows();
	int eta_size = (*conservative).cols();
	for (int i = 0; i < xi_size; i++)
	{
		for (int j = 0; j < eta_size; j++)
		{
			double dt = maxCFL * calcTimeStep(i, j, conservative); //local time step
			Dxi = (*xi_fluxes)(i+1,j) * grid->getSxi(i + 1, j) - (*xi_fluxes)(i,j) * grid->getSxi(i, j);
			Deta = (*eta_fluxes)(i,j+1) * grid->getSeta(i, j + 1) - (*eta_fluxes)(i,j) * grid->getSeta(i, j);
			(*conservative)(i,j) = (*conservative)(i,j) - dt / grid->getVolume(i, j) * (Dxi + Deta);
		}
	}
}