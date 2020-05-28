#pragma once
#include "TimeStepper.h"
class TimeStepperLocal :
	public TimeStepper
{
public:
	TimeStepperLocal();
	~TimeStepperLocal();

	void execute(StateMatrix2D * conservative, StateMatrix2D *xi_fluxes, StateMatrix2D *eta_fluxes);
};

