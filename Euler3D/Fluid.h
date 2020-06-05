#pragma once
#include "GlobalTypes.h"

#include <Eigen/Dense>

class Fluid
{
public:
	Fluid();
	virtual ~Fluid();

	virtual void setGamma(double gamma) = 0;
	virtual void setCp(double cp) = 0;

	double calcKineticState(const StateVector & primcons);
	double calcKineticU(const DirVector & u);

	virtual double getSoundSpeed(const StateVector &prim) = 0;
	virtual double getPressure(const StateVector &prim) = 0;
	virtual double getTemperature(const StateVector &prim) = 0;
	virtual double getEnthalpy(const StateVector &prim) = 0;
	virtual double getGamma() = 0;
	virtual double getCp() = 0;

	virtual void setAllFromGammaCp() = 0;

	//Utilities for conversion between variable types
	virtual StateVector prim2cons(const StateVector &p) = 0;
	virtual StateVector cons2prim(const StateVector &c) = 0;
	virtual StateVector user2cons(double p, const DirVector &uvec, double T) = 0;
	virtual double calcPcons(const StateVector & c) = 0;
	virtual double calcMacons(const StateVector & c) = 0;
	virtual double calcSoundSpeedcons(const StateVector & c) = 0;
	virtual double calcPprim(const StateVector & prim) = 0;
	virtual StateVector calcPhysFlux(const StateVector & c, DirVector &n) = 0;
};

