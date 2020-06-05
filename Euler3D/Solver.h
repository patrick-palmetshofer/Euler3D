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
	std::array<std::unique_ptr<Flux>, 3> fluxes;
	std::vector<std::unique_ptr<Boundary>> boundaries;

	Eigen::Array<StateTensor *, 3, 1> flux_tensors;

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
	StateVector cons_inlet;
	StateVector cons_initial;
	//StateTensor primitive;

	//StateTensor xi_face_prim;
	//StateTensor eta_face_prim;

	StateTensor conservative;
	StateTensor old_conservative;

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
	std::vector<std::unique_ptr<Boundary>> * setBoundaries(std::vector<std::unique_ptr<Boundary>>& new_boundaries);
	Eigen::Array<std::unique_ptr<Flux>,3, 1> * setFluxes(Eigen::Array<std::unique_ptr<Flux>,3,1> &new_fluxes);

	void crossPopulatePointers();

	void setTimeWriteInterval(double t) { time_write_interval = t; };
	void setIterWriteInterval(int i) { iter_write_interval = i; };

	//Utilities for Sod problem tests
	void setSodXInitial();
	void setSodYInitial();

	//Set initial and boundary condition values for initialization
	void setConsInlet(double p, const DirVector &uvec, double T);
	void setConsInitial(double p, const DirVector &uvec, double T);

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
	//StateVector calcFlux();
	//StateVector spaceDisc(int i, int j, int k);
	//StateVector spaceDisc(const StateVector & c, DirVector &n);
	//StateVector calcFlux(const StateVector & c, DirVector &n);
	//StateVector calcPhysFluxesXi();
	//StateVector mainloop();
	//StateVector calcDissip();
	//void physicalFluxRight(const StateVector &c);
	//void physicalFluxUp(const StateVector &c);
	//void RoeDissipRight(int i, int j, int k);
	//StateVector calcPhysFlux(const StateVector & prim_left, const StateVector & prim_right, DirVector &n);
};

