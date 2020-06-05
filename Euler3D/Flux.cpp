#include "Flux.h"

using namespace Euler;

Flux::Flux(int new_dim) : dim(new_dim)
{
	dim_arr.fill(0);
	dim_arr[dim] = 1;
}


Flux::~Flux()
{
}

void Flux::setConservative(StateTensor * cons)
{
	conservative = cons;
	max_inds = grid->getnComponentCells() + dim_arr;
	Euler::resize(fluxes, max_inds);
	if (max_inds[dim] < 3)
		max_inds[dim] = 0;
}

void Flux::calcFluxes()
{
	IndArray max_temp = max_inds - dim_arr;
	for (int i = dim_arr[0]; i < max_temp[0]; i++)
	{
		for (int j = dim_arr[1]; j < max_temp[1]; j++)
		{
			for (int k = dim_arr[2]; k < max_temp[2]; k++)
			{
				leftrightstates = reconstruct->reconstructStates(i, j, k, dim_arr);
				DirVector n0 = grid->getnVec(i, j, k, dim);
				DirVector n = Euler::swap(n0, 0, dim);

				leftrightstates = { swap(leftrightstates.first), swap(leftrightstates.second) };
				flux = calcFlux(leftrightstates, n);

				if (Euler::checkNaN(&flux))
					throw;

				StateVector fl = swap(flux);
				fluxes(i)(j)(k) = fl;
			}
		}
	}
}

StateVector Flux::calcFlux(std::pair<StateVector, StateVector> leftrightstates, DirVector &n)
{
	StateVector leftflux = fluid->calcPhysFlux(leftrightstates.first, n);
	StateVector rightflux = fluid->calcPhysFlux(leftrightstates.second, n);
	StateVector dissip = calcDissip(leftrightstates, n);
	return 0.5*(leftflux + rightflux - dissip);
}
