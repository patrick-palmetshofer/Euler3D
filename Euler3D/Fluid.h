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

	virtual double getSoundSpeed(const StateVector2D &prim) = 0;
	virtual double getPressure(const StateVector2D &prim) = 0;
	virtual double getTemperature(const StateVector2D &prim) = 0;
	virtual double getEnthalpy(const StateVector2D &prim) = 0;
	virtual double getGamma() = 0;
	virtual double getCp() = 0;

	virtual void setAllFromGammaCp() = 0;

	//Utilities for conversion between variable types
	virtual StateVector2D prim2cons(const StateVector2D &p) = 0;
	virtual StateVector2D cons2prim(const StateVector2D &c) = 0;
	virtual StateVector2D user2cons(double p, double u, double v, double T) = 0;
	virtual double calcPcons(const StateVector2D & c) = 0;
	virtual double calcMacons(const StateVector2D & c) = 0;
	virtual double calcSoundSpeedcons(const StateVector2D & c) = 0;
	virtual double calcPprim(const StateVector2D & prim) = 0;
	virtual StateVector2D calcPhysFlux(const StateVector2D & c, Eigen::Vector2d& n) = 0;
	virtual StateVector2D calcPhysFlux(const StateVector2D & c, double nx, double ny) = 0;
};

