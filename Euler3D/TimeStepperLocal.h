#pragma once
#include "TimeStepper.h"
class TimeStepperLocal :
	public TimeStepper
{
public:
	TimeStepperLocal();
	~TimeStepperLocal();

	void execute(StateTensor * conservative, Eigen::Array<StateTensor *, 3, 1>);
};

