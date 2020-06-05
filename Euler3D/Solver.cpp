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
	Euler::resize(conservative, grid->getnComponentCells());
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

Eigen::Array<std::unique_ptr<Flux>,3, 1>* Solver::setFluxes(Eigen::Array<std::unique_ptr<Flux>,3,1>& new_fluxes)
{
	for (int dim = 0; dim < 3; ++dim)
	{
		fluxes[dim] = std::move(new_fluxes[dim]);
	}
	return nullptr;
}

void Solver::crossPopulatePointers()
{
	stepper->setGrid(grid.get());
	stepper->setFluid(fluid.get());
	stepper->setCFL(maxCFL);

	reconstruct->setGrid(grid.get());

	solution->setFluid(fluid.get());

	for (int dim = 0; dim < 3; ++dim)
	{
		fluxes[dim]->setFluid(fluid.get());
		fluxes[dim]->setGrid(grid.get());
		fluxes[dim]->setReconstruct(reconstruct.get());
	}

	for (auto &b : boundaries)
		b->setGrid(grid.get());
}

//Sets conservative variables at inlet, takes user primitive variables
void Solver::setConsInlet(double p, const DirVector &uvec, double T)
{
	cons_inlet = fluid->user2cons(p, uvec, T);
}

//Sets conservative variables at initial time, takes user primitive variables
//Sets all cells in domain to the specified value
void Solver::setConsInitial(double p, const DirVector &uvec, double T)
{
	p_infty = p;
	cons_initial = fluid->user2cons(p, uvec, T);
	setInitialCondition();
}

//Takes cons-initial Solver parameter and sets all cells in domain to these conditions
void Solver::setInitialCondition()
{
	Euler::fill(conservative, cons_initial);
}

double Solver::getTime()
{
	return time;
}

void Solver::allocateConservative()
{
	Euler::resize(conservative, grid->getnComponentCells());
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
	for (int dim = 0; dim < 3; ++dim)
	{
		if (!fluxes[dim])
			throw;
		flux_tensors[dim] = fluxes[dim]->get();
	}
	

	allocateConservative();
	setInitialCondition();

	for (auto &b : boundaries)
	{
		for (int dim = 0; dim < 3; ++dim)
		{
			if (b->getDim() == fluxes[dim]->getDim())
				b->setFlux(fluxes[dim].get());
		}
		b->setConservative(&conservative);
		b->setFluid(fluid.get());
		b->init();
	}
	reconstruct->setConservative(&conservative);
	reconstruct->setGrid(grid.get());
	for (int dim = 0; dim < 3; ++dim)
	{
		fluxes[dim]->setConservative(&conservative);
	}
	for (int i = 0; i <= maxIter; i++)
	{
		//StateVector calcFluxMUSCL(const StateVector & c_left_left, const StateVector & c_left, const StateVector & c_right, const StateVector & c_right_right, DirVector &n);
		//StateVector calcFluxMUSCL(const StateVector & c_left_left, const StateVector & c_left, const StateVector & c_right, const StateVector & c_right_right, double Sx_left_left, double Sx_left, double Sx_right, double Sx_right_right, DirVector &n);
		for (auto &b : boundaries)
			b->apply();

		//Calculate fluxes
		//#pragma omp parallel for
		for (int dim = 0; dim < 3; ++dim)
		{
			fluxes[dim]->calcFluxes();
			if (Euler::checkNaN(fluxes[dim]->get()))
				throw;
		}

		old_conservative = conservative;

		stepper->execute(&conservative, flux_tensors);

		if (Euler::checkNaN(&conservative))
		{
			solution->writeSolution(old_conservative, i-1);
			throw;
		}

		solution->calcResidualsL2(old_conservative,conservative);
		solution->calcResidualsLinfty(old_conservative,conservative);
		
		if (i % iter_write_interval == 0 || Euler::checkNaN(&conservative))
		{
			//solution->writeSolution(conservative, i);
			solution->writeSolution("solution/newbump.cgns", grid.get(), conservative);
		}
	}
}