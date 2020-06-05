#include "TimeStepperLocal.h"



TimeStepperLocal::TimeStepperLocal()
{
}


TimeStepperLocal::~TimeStepperLocal()
{
}

void TimeStepperLocal::execute(StateTensor * conservative, Eigen::Array<StateTensor *,3,1> fluxes)
{
	for (int i = 0; i < conservative->size(); i++)
	{
		for (int j = 0; j < (*conservative)(0).size(); j++)
		{
			for (int k = 0; k < (*conservative)(0)(0).size(); k++)
			{
				double dt = maxCFL * calcTimeStep(i, j, k, conservative); //local time step
				IndArray ind(i, j, k);
				diffs.fill(0);
				for (int dim = 0; dim < ind.size(); ++dim)
				{
					IndArray indplus = ind;
					indplus[dim] = ind[dim]+1;

					DirVector n = grid->getnVec(ind, dim);
					DirVector nplus = grid->getnVec(indplus, dim);

					StateVector fluxplus = (*(fluxes(dim)))(indplus[0])(indplus[1])(indplus[2]);
					StateVector fluxzero = (*(fluxes(dim)))(ind[0])(ind[1])(ind[2]);

					double fplus = grid->getFaceArea(ind, dim);
					double fzero = grid->getFaceArea(indplus, dim);

					diffs += fluxplus * fplus - fluxzero * fzero;
				}

				double vol = grid->getVolume(i,j,k);
				(*conservative)(i)(j)(k) -= dt / vol * diffs;
			}
		}
	}
}