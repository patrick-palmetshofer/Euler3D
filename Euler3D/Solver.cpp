#include "Solver.h"
#include <fstream>
#include <iostream>

//Ugly, hard-coded initial conditions for Sod problem

//Constructor for solver class. Initializes time and MUSCL parameters
Solver::Solver() : Solver(1,1.0/3.0)
{
}

Solver::Solver(double eps, double kappa)
{
	time = 0;
	limit = true;
	//reconstruct_eps = eps;
	//reconstruct_kappa = kappa;// 1.0 / 3.0;// 1.0 / 3.0;// 1.0 / 3.0;// 1.0 / 3.0;// 1 / 3;
}

Solver::Solver(std::string filename, double eps, double kappa) : Solver(eps, kappa)
{
	grid->readGridPro(filename);
	initSizeFromGrid();
}

//Consolidated constructor and file reader
Solver::Solver(std::string filename) : Solver()
{
	grid->readGridPro(filename);
	initSizeFromGrid();
}

void Solver::initSizeFromGrid()
{
	conservative.resize(grid->getnxiCells(), grid->getnetaCells());
}

//Destructor. As everything is implemented using STL containers, data deallocation is handled by STL
Solver::~Solver()
{
}

std::vector<std::unique_ptr<Boundary>>* Solver::setBoundaries(std::vector<std::unique_ptr<Boundary>>& new_boundaries)
{
	boundaries.clear(); 
	boundaries.reserve(new_boundaries.size());
	for (auto &b : new_boundaries) 
	{ 
		boundaries.push_back(std::move(b));
	}
	return &boundaries;
}

void Solver::crossPopulatePointers()
{
	stepper->setGrid(grid.get());
	stepper->setFluid(fluid.get());
	stepper->setCFL(maxCFL);

	reconstruct->setGrid(grid.get());

	solution->setFluid(fluid.get());

	xi_fluxes->setFluid(fluid.get());
	xi_fluxes->setGrid(grid.get());
	xi_fluxes->setReconstruct(reconstruct.get());

	eta_fluxes->setFluid(fluid.get());
	eta_fluxes->setGrid(grid.get());
	eta_fluxes->setReconstruct(reconstruct.get());

	for (auto &b : boundaries)
		b->setGrid(grid.get());
}

//Sets conservative variables at inlet, takes user primitive variables
void Solver::setConsInlet(double p, double u, double v, double T)
{
	cons_inlet = fluid->user2cons(p, u, v, T);
}

//Sets conservative variables at initial time, takes user primitive variables
//Sets all cells in domain to the specified value
void Solver::setConsInitial(double p, double u, double v, double T)
{
	p_infty = p;
	cons_initial = fluid->user2cons(p, u, v, T);
	setInitialCondition();
}

//Takes cons-initial Solver parameter and sets all cells in domain to these conditions
void Solver::setInitialCondition()
{
	conservative.fill(cons_initial);
}

double Solver::getTime()
{
	return time;
}

void Solver::allocateConservative()
{
	conservative.resize(grid->getnxiCells(), grid->getnetaCells());
}

void Solver::solve()
{
	if (!grid)
		throw;
	if (!stepper)
		throw;
	if (!reconstruct)
		throw;
	if (!fluid)
		throw;
	if (!solution)
		throw;
	if (!xi_fluxes)
		throw;
	if (!eta_fluxes)
		throw;

	allocateConservative();
	setInitialCondition();

	for (auto &b : boundaries)
	{
		if (b->getDim() == xi_fluxes->getDim())
			b->setFlux(xi_fluxes.get());
		else if (b->getDim() == eta_fluxes->getDim())
			b->setFlux(eta_fluxes.get());
		else
			throw;
		b->setConservative(&conservative);
		b->setFluid(fluid.get());
		b->init();
	}
	reconstruct->setConservative(&conservative);
	reconstruct->setGrid(grid.get());
	xi_fluxes->setConservative(&conservative);
	eta_fluxes->setConservative(&conservative);
	for (int i = 0; i <= maxIter; i++)
	{
		//StateVector2D calcFluxMUSCL(const StateVector2D & c_left_left, const StateVector2D & c_left, const StateVector2D & c_right, const StateVector2D & c_right_right, double nx, double ny);
		//StateVector2D calcFluxMUSCL(const StateVector2D & c_left_left, const StateVector2D & c_left, const StateVector2D & c_right, const StateVector2D & c_right_right, double Sx_left_left, double Sx_left, double Sx_right, double Sx_right_right, double nx, double ny);
		for (auto &b : boundaries)
			b->apply();

		//Calculate fluxes
		xi_fluxes->calcFluxes();
		eta_fluxes->calcFluxes();

		//Euler::checkNaN(xi_fluxes->get());
		//Euler::checkNaN(eta_fluxes->get());

		old_conservative = conservative;

		stepper->execute(&conservative, xi_fluxes->get(), eta_fluxes->get());

		if (Euler::checkNaN(&conservative))
		{
			solution->writeSolution(old_conservative, i-1);
			throw;
		}

		solution->calcResidualsL2(old_conservative,conservative);
		solution->calcResidualsLinfty(old_conservative,conservative);
		
		if (i % iter_write_interval == 0 || Euler::checkNaN(&conservative))
			solution->writeSolution(conservative,i);
	}
}