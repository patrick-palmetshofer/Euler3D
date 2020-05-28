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

double PerfectGas::getSoundSpeed(const StateVector2D &prim)
{
	double T = getTemperature(prim);
	return std::sqrt(gamma*R*T);
}

double PerfectGas::getPressure(const StateVector2D &prim)
{
	return prim[0]*R*getTemperature(prim);
}

double PerfectGas::getTemperature(const StateVector2D &prim)
{
	return prim[3]/cv;
}

double PerfectGas::getEnthalpy(const StateVector2D &prim)
{
	return gamma*prim[3];
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
StateVector2D PerfectGas::prim2cons(const StateVector2D & p)
{
	StateVector2D c;
	c[0] = p[0];
	c[1] = p[0] * p[1];
	c[2] = p[0] * p[2];
	c[3] = p[0] * p[3];
	return c;
}

StateVector2D PerfectGas::cons2prim(const StateVector2D & c)
{
	StateVector2D p;
	p[0] = c[0];
	p[1] = c[1] / c[0];
	p[2] = c[2] / c[0];
	p[3] = c[3] / c[0];
	return p;
}

StateVector2D PerfectGas::user2cons(double p, double u, double v, double T)
{
	//Easy input variables to conservative
	double e = cp / gamma * T;
	double et = e + 0.5*(u*u + v * v);
	double rho = p / ((gamma - 1)*e);
	StateVector2D prim = { rho,u,v,et };
	StateVector2D cons = prim2cons(prim);
	return cons;
}

//Functions to calculate physical properties (and charactersitic variables) from conservative/primite variables
double PerfectGas::calcPcons(const StateVector2D & c)
{
	double p = (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0]);
	return p;
}

double PerfectGas::calcMacons(const StateVector2D & c)
{
	double Ma = std::sqrt(c[1] * c[1] + c[2] * c[2]) / (c[0] * calcSoundSpeedcons(c));
	return Ma;
}

double PerfectGas::calcSoundSpeedcons(const StateVector2D &c)
{
	double sound = std::sqrt((gamma - 1)*gamma*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0]) / c[0]);
	return sound;
}

double PerfectGas::calcPprim(const StateVector2D & prim)
{
	double p = prim[0] * (gamma - 1)*(prim[3] - 0.5*(prim[1] * prim[1] + prim[2] * prim[2]));
	return p;
}

//Calculate physical flux without numerical dissipation scheme
StateVector2D PerfectGas::calcPhysFlux(const StateVector2D & c, Eigen::Vector2d& n)
{
	double &nx = n[0];
	double &ny = n[1];
	return calcPhysFlux(c, nx, ny);
}

StateVector2D PerfectGas::calcPhysFlux(const StateVector2D & c, double nx, double ny)
{
	StateVector2D flux;
	flux[0] = c[1] * nx + c[2] * ny;
	flux[1] = c[1] * c[1] / c[0] * nx + c[1] * c[2] / c[0] * ny + (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0])*nx;
	flux[2] = c[2] * c[1] / c[0] * nx + c[2] * c[2] / c[0] * ny + (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0])*ny;
	flux[3] = (gamma*c[3] - 0.5*(gamma - 1)*(c[1] * c[1] + c[2] * c[2]) / c[0])*(c[1] * nx + c[2] * ny) / c[0];
	return flux;
}