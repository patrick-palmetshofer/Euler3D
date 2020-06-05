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

	double getSoundSpeed(const StateVector &prim);
	double getPressure(const StateVector &prim);
	double getTemperature(const StateVector &prim);
	double getEnthalpy(const StateVector &prim);
	double getGamma();
	double getCp();

	void setAllFromGammaCp();

	//Utilities for conversion between variable types
	StateVector prim2cons(const StateVector &p);
	StateVector cons2prim(const StateVector &c);
	StateVector user2cons(double p, const DirVector &u, double T);
	double calcPcons(const StateVector & c);
	double calcMacons(const StateVector & c);
	double calcSoundSpeedcons(const StateVector & c);
	double calcPprim(const StateVector & prim);
	StateVector calcPhysFlux(const StateVector & c, DirVector &n);
};
