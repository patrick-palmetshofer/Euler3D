#pragma once
#include "GlobalTypes.h"
#include "Grid.h"
#include "Flux.h"
class TimeStepper
{
protected:
	Grid * grid;
	Fluid * fluid;

	//Calculate time step for global time stepping
	double calcTimeStep(StateMatrix2D * conservative);
	//Method calls:
	//Calculate timestep for local timestepping
	double calcTimeStep(int i, int j, StateMatrix2D * conservative);

	StateVector2D p;
	double maxCFL = 0.5;

	StateVector2D Dxi, Deta;

public:
	TimeStepper();
	virtual ~TimeStepper();

	inline void setFluid(Fluid * new_fluid) { fluid = new_fluid; };
	inline void setGrid(Grid * new_grid) { grid = new_grid; };
	void setCFL(double cfl) { maxCFL = cfl; };

	virtual void execute(StateMatrix2D * conservative, StateMatrix2D * xi_fluxes, StateMatrix2D * eta_fluxes) = 0;
};

