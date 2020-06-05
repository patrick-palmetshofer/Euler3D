#include "PerfectGas.h"
#include <cmath>


//W in g/mol; cp in J/kgK
PerfectGas::PerfectGas(double W,double cp)
{
	this->W = W;
	this->cp = cp;

	R = 1e3*Rm / W;
	cv = cp - R;
	gamma = cp / cv;
}


void PerfectGas::setAllFromGammaCp()
{
	cv = cp / gamma;
	R = cp - cv;
	W = 1e3*Rm / R;
}

PerfectGas::~PerfectGas()
{
}

double PerfectGas::getSoundSpeed(const StateVector &prim)
{
	double T = getTemperature(prim);
	return std::sqrt(gamma*R*T);
}

double PerfectGas::getPressure(const StateVector &prim)
{
	return prim[4]*R*getTemperature(prim);
}

double PerfectGas::getTemperature(const StateVector &prim)
{
	return prim[4]/cv;
}

double PerfectGas::getEnthalpy(const StateVector &prim)
{
	return gamma*prim[4];
}

double PerfectGas::getGamma()
{
	return gamma;
}

double PerfectGas::getCp()
{
	return cp;
}

//Functions to convert between conservative and primitive variables
StateVector PerfectGas::prim2cons(const StateVector & p)
{
	StateVector c;
	c[0] = p[0];
	for (int i = 1; i < 5; ++i)
		c[i] = p[0] * p[i];
	return c;
}

StateVector PerfectGas::cons2prim(const StateVector & c)
{
	StateVector p;
	p[0] = c[0];
	for (int i = 1; i < 5; ++i)
		p[i] = c[i] / c[0];
	return p;
}

StateVector PerfectGas::user2cons(double p, const DirVector &u, double T)
{
	//Easy input variables to conservative
	double e = cp / gamma * T;
	double et = e + calcKineticU(u);
	double rho = p / ((gamma - 1)*e);
	StateVector prim;
	prim << rho, u[0], u[1], u[2], et;
	StateVector cons = prim2cons(prim);
	return cons;
}

//Functions to calculate physical properties (and charactersitic variables) from conservative/primite variables
double PerfectGas::calcPcons(const StateVector & c)
{
	double p = (gamma - 1)*(c[4] - calcKineticState(c) / c[0]);
	return p;
}

double PerfectGas::calcMacons(const StateVector & c)
{
	double Ma = std::sqrt(calcKineticState(c)) / (c[0] * calcSoundSpeedcons(c));
	return Ma;
}

double PerfectGas::calcSoundSpeedcons(const StateVector &c)
{
	double sound = std::sqrt((gamma - 1)*gamma*(c[4] - calcKineticState(c) / c[0]) / c[0]);
	return sound;
}

double PerfectGas::calcPprim(const StateVector & prim)
{
	double p = prim[0] * (gamma - 1)*(prim[4] - calcKineticState(prim));
	return p;
}

//Calculate physical flux without numerical dissipation scheme
StateVector PerfectGas::calcPhysFlux(const StateVector & c, DirVector &n)
{
	StateVector flux;
	flux[0] = c.segment(1,3).matrix().dot(n);
	flux.segment(1,3) = c.segment(1,3) / c[0] * flux[0] + (gamma - 1)*(c[4] - calcKineticState(c) / c[0])*n.array();
	flux[4] = (gamma*c[4] - 0.5*(gamma - 1)*(calcKineticState(c)) / c[0])*flux[0] / c[0];
	return flux;

	//flux[1] = c[1] / c[0] * flux[0] + (gamma - 1)*(c[4] - 0.5*(calcKineticState(c)) / c[0])*n(0);
	//flux[2] = c[2] / c[0] * flux[0] + (gamma - 1)*(c[4] - 0.5*(calcKineticState(c)) / c[0])*n(1);
	//flux[3] = c[3] / c[0] * flux[0] + (gamma - 1)*(c[4] - 0.5*(calcKineticState(c)) / c[0])*n(2);
}