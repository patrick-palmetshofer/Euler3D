#include "Fluid.h"



Fluid::Fluid()
{
}


Fluid::~Fluid()
{
}

double Fluid::calcKineticState(const StateVector & primcons)
{
	return calcKineticU(primcons.segment(1,3).matrix());
}

double Fluid::calcKineticU(const DirVector & u)
{
	return 0.5*u.squaredNorm();
}
