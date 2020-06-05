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

	void execute(StateTensor * conservative, Eigen::Array<StateTensor *, 3, 1> fluxes);
};

