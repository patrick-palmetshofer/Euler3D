#pragma once
#include "GlobalTypes.h"
#include "Fluid.h"
#include <Eigen/Dense>

class PerfectGas : 
	public Fluid
{
private:
	const double Rm = 8.31446261815324;
	//Perfect gas constants
	double gamma = 1.4;
	double cp = 1005;
	double W;
	double cv;
	double R;

public:
	PerfectGas(double gamma, double cp);
	~PerfectGas();

	void setGamma(double gamma) { this->gamma = gamma; };
	void setCp(double cp) { this->cp = cp; };

	double getSoundSpeed(const StateVector2D &prim);
	double getPressure(const StateVector2D &prim);
	double getTemperature(const StateVector2D &prim);
	double getEnthalpy(const StateVector2D &prim);
	double getGamma();
	double getCp();

	void setAllFromGammaCp();

	//Utilities for conversion between variable types
	StateVector2D prim2cons(const StateVector2D &p);
	StateVector2D cons2prim(const StateVector2D &c);
	StateVector2D user2cons(double p, double u, double v, double T);
	double calcPcons(const StateVector2D & c);
	double calcMacons(const StateVector2D & c);
	double calcSoundSpeedcons(const StateVector2D & c);
	double calcPprim(const StateVector2D & prim);
	StateVector2D calcPhysFlux(const StateVector2D & c, Eigen::Vector2d& n);
	StateVector2D calcPhysFlux(const StateVector2D & c, double nx, double ny);
};
