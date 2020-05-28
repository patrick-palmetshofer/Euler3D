#pragma once
#include "GlobalTypes.h"
#include "Fluid.h"
#include "Reconstruct.h"

class Flux
{
private:
	std::pair<StateVector2D, StateVector2D> leftrightstates;
	StateVector2D flux;
protected:
	Reconstruct * reconstruct;
	Fluid * fluid;
	Grid * grid;

	//Dimension of the flux object: 1 for Xi, 2 for Eta
	const int dim;

	//Fluxes at Faces and Conservative variables of all cells
	StateMatrix2D fluxes;
	StateMatrix2D *conservative;

	template<typename T>
	T swap(T &data)
	{
		return Euler::swap(data, dim+1);
	}

	//std::pair<StateVector2D, StateVector2D> reconstruct(int i, int j);

public:
	Flux(int new_dim);
	virtual ~Flux();

	void setFluid(Fluid * new_fluid) { fluid = new_fluid; };
	void setReconstruct(Reconstruct * new_reconstruct) { reconstruct = new_reconstruct; };
	void setGrid(Grid * new_grid) { grid = new_grid; };
	void setConservative(StateMatrix2D * cons);
	void setBoundaryFlux(int i, int j, StateVector2D flux) { fluxes(i,j) = flux; };

	StateMatrix2D * get() { return &fluxes; };

	//Calculate fluxes of all faces.
	//Xi=const faces

	void calcFluxes();


	int getDim() { return dim; };

	StateVector2D calcFlux(std::pair<StateVector2D, StateVector2D> leftrightstates, double nx, double ny);

	//Flux calculation
	//Previous method calls:
	//Calculate physical fluxes without numerical dissipation
	StateVector2D calcPhysFlux(const StateVector2D & c, double nx, double ny);
	//Calculate numerical dissipation
	virtual StateVector2D calcDissip(std::pair<StateVector2D, StateVector2D> leftrightstates, double nx, double ny) = 0;
};

