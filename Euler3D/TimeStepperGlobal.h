#pragma once
#include "TimeStepper.h"

class TimeStepperGlobal :
	public TimeStepper
{
private:
	double dt;
	double time;
public:
	TimeStepperGlobal();
	~TimeStepperGlobal();

	void execute(StateMatrix2D * conservative, StateMatrix2D * xi_fluxes, StateMatrix2D * eta_fluxes);
};

