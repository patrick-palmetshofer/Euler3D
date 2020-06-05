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
	double calcTimeStep(StateTensor * conservative);
	//Method calls:
	//Calculate timestep for local timestepping
	double calcTimeStep(int i, int j, int k, StateTensor * conservative);

	StateVector p;
	double maxCFL = 0.5;

	StateVector diffs;

public:
	TimeStepper();
	virtual ~TimeStepper();

	inline void setFluid(Fluid * new_fluid) { fluid = new_fluid; };
	inline void setGrid(Grid * new_grid) { grid = new_grid; };
	void setCFL(double cfl) { maxCFL = cfl; };

	virtual void execute(StateTensor * conservative, Eigen::Array<StateTensor *, 3, 1>) = 0;
};

