#pragma once
#include "GlobalTypes.h"
#include "Grid.h"
#include "TimeStepper.h"
#include "Fluid.h"
#include "Solution.h"
#include "Boundary.h"

#include <memory>

class Solver
{
private:
	std::unique_ptr<Grid> grid;
	std::unique_ptr<TimeStepper> stepper;
	std::unique_ptr<Reconstruct> reconstruct;
	std::unique_ptr<Fluid> fluid;
	std::unique_ptr<Solution> solution;
	std::unique_ptr<Flux> xi_fluxes;
	std::unique_ptr<Flux> eta_fluxes;
	std::vector<std::unique_ptr<Boundary>> boundaries;

	double time_write_interval;
	int iter_write_interval;

	//Use limiter?
	bool limit;

	double implicit;

	//Maximum CFL number, use 0.5 for sure stability
	double maxCFL = 0.1;
	double p_infty;

	//Time step and current time in simulation
	double dt;
	double time;

	double maxTime;
	double maxIter;

	//conservative variables for Inlet and initial condition
	StateVector2D cons_inlet;
	StateVector2D cons_initial;
	//StateMatrix2D primitive;

	//StateMatrix2D xi_face_prim;
	//StateMatrix2D eta_face_prim;

	StateMatrix2D conservative;
	StateMatrix2D old_conservative;

	//Treatment of Boundary conditions. Hard coded :(
	//First-order (constant boundary) conditions
	void setBoundaryInletLeft();
	void setBoundaryLowerWall();
	void setBoundaryUpperLowerWalls();
	void setWalls();
	void setBoundaryOutlet();
	void setBoundaryUpperOutlet();
	//LODI Boundary conditions
	void setCharacteristicBoundaryRightOutlet();
	void setCharacteristicBoundaryUpperOutlet();

public:

	void setImplicit(double imp) { implicit = imp; };
	void setMaxTime(double t) { maxTime = t; };
	void setMaxIter(double iter) { maxIter = iter; };
	void setCFL(double c) { maxCFL = c; };

	template <typename T> T * set(std::unique_ptr<T>& source, std::unique_ptr<T>& target) { target = std::move(source); return target.get(); };

	TimeStepper * setTimeStepper(std::unique_ptr<TimeStepper>& new_stepper) { return set(new_stepper, stepper); };
	Reconstruct * setReconstruct(std::unique_ptr<Reconstruct>& new_reconstruct) { return set(new_reconstruct, reconstruct); };
	Grid * setGrid(std::unique_ptr<Grid>& new_grid) { return set(new_grid, grid); };
	Fluid * setFluid(std::unique_ptr<Fluid>& new_fluid) { return set(new_fluid, fluid); };
	Solution * setSolution(std::unique_ptr<Solution>& new_solution) { return set(new_solution, solution); };
	Flux * setXiFluxes(std::unique_ptr<Flux>& new_xi_fluxes) { return set(new_xi_fluxes, xi_fluxes); };
	Flux * setEtaFluxes(std::unique_ptr<Flux>& new_eta_fluxes) { return set(new_eta_fluxes, eta_fluxes); };
	std::vector<std::unique_ptr<Boundary>> * setBoundaries(std::vector<std::unique_ptr<Boundary>>& new_boundaries);

	void crossPopulatePointers();

	void setTimeWriteInterval(double t) { time_write_interval = t; };
	void setIterWriteInterval(int i) { iter_write_interval = i; };

	//Utilities for Sod problem tests
	void setSodXInitial();
	void setSodYInitial();

	//Set initial and boundary condition values for initialization
	void setConsInlet(double p, double u, double v, double T);
	void setConsInitial(double p, double u, double v, double T);

	//Fill all cells with initial condition
	void setInitialCondition();

	//Constructors/destructors
	Solver();
	Solver(double eps, double kappa);
	Solver(std::string filename, double eps, double kappa);
	Solver(std::string filename);
	void initSizeFromGrid();
	~Solver();

	//Get current time for global time stepping. Returns 0 for local time steps
	double getTime();

	void allocateConservative();

	//limiter control
	void enableLimiter() {
		limit = true;
	};
	void disableLimiter() {
		limit = false;
	};

	void solve();

	// Deprecated prototypes
	//StateVector2D calcFlux();
	//StateVector2D spaceDisc(int i, int j);
	//StateVector2D spaceDisc(const StateVector2D & c, double nx, double ny);
	//StateVector2D calcFlux(const StateVector2D & c, double nx, double ny);
	//StateVector2D calcPhysFluxesXi();
	//StateVector2D mainloop();
	//StateVector2D calcDissip();
	//void physicalFluxRight(const StateVector2D &c);
	//void physicalFluxUp(const StateVector2D &c);
	//void RoeDissipRight(int i, int j);
	//StateVector2D calcPhysFlux(const StateVector2D & prim_left, const StateVector2D & prim_right, double nx, double ny);
};

