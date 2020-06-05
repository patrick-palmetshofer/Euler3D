#pragma once
#include "GlobalTypes.h"
#include "Fluid.h"
#include "Reconstruct.h"

class Flux
{
private:
	std::pair<StateVector, StateVector> leftrightstates;
	StateVector flux;
protected:
	Reconstruct * reconstruct;
	Fluid * fluid;
	Grid * grid;

	//Dimension of the flux object: 1 for Xi, 2 for Eta
	const int dim;

	Eigen::Array3i dim_arr;
	Eigen::Array3i max_inds;

	//Fluxes at Faces and Conservative variables of all cells
	StateTensor fluxes;
	StateTensor *conservative;

	template<typename T>
	T swap(T &data)
	{
		return Euler::swap(data, dim+1);
	}

	//std::pair<StateVector, StateVector> reconstruct(int i, int j, int k);

public:
	Flux(int new_dim);
	virtual ~Flux();

	void setFluid(Fluid * new_fluid) { fluid = new_fluid; };
	void setReconstruct(Reconstruct * new_reconstruct) { reconstruct = new_reconstruct; };
	void setGrid(Grid * new_grid) { grid = new_grid; };
	void setConservative(StateTensor * cons);
	void setBoundaryFlux(int i, int j, int k, StateVector flux) { fluxes(i)(j)(k) = flux; };
	void setBoundaryFlux(IndArray &ind, StateVector flux) { fluxes(ind[0])(ind[1])(ind[2]) = flux; };

	StateTensor * get() { return &fluxes; };

	//Calculate fluxes of all faces.
	//Xi=const faces

	void calcFluxes();


	int getDim() { return dim; };

	StateVector calcFlux(std::pair<StateVector, StateVector> leftrightstates, DirVector &n);

	//Flux calculation
	//Previous method calls:
	//Calculate physical fluxes without numerical dissipation
	StateVector calcPhysFlux(const StateVector & c, DirVector &n);
	//Calculate numerical dissipation
	virtual StateVector calcDissip(std::pair<StateVector, StateVector> leftrightstates, DirVector &n) = 0;
};

