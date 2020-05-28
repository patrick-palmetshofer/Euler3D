#pragma once
#include "GlobalTypes.h"
#include "Grid.h"
#include "Fluid.h"
#include "Flux.h"

class Boundary
{
protected:
	Grid *grid;
	Fluid *fluid;
	StateMatrix2D * conservative;
	Flux * flux;

	int dim;
	int dir;
	std::vector<std::array<int, 2>> cell_inds;

	StateVector2D ghost_state;

	template<typename T>
	T swap(T &data)
	{
		return Euler::swap(data, dim);
	}
public:
	Boundary();
	virtual ~Boundary();

	void setBottom();
	void setTop();
	void setRight();
	void setLeft();

	void init();

	void setGrid(Grid * new_grid) { grid = new_grid; };
	void setFluid(Fluid * newfluid) { fluid = newfluid; };
	void setFlux(Flux * newflux) { flux = newflux; };
	void setConservative(StateMatrix2D * cons) { conservative = cons; };

	int getDim() { return dim; };
	int getDir() { return dir; };

	virtual StateVector2D getGhostState(std::array<int,2> ind) = 0;

	void apply();
};

class BoundaryLODI :
	public Boundary
{
protected:
	double L1 = 0, L2 = 0, L3 = 0, L5 = 0;
public:
	virtual void calcWaveStrengths(double u, double sound, double rho, double dx, double drho, double dp, double du, StateVector2D & c, StateVector2D & cneg) = 0;
	StateVector2D getGhostState(std::array<int, 2> ind);
};

class SupersonicOutlet :
	public Boundary
{
public:
	StateVector2D getGhostState(std::array<int, 2> ind);
};

class Outlet :
	public BoundaryLODI
{
private:
	double p_infty;
public:
	Outlet(double new_pinfty) : p_infty(new_pinfty) {};
	void setPinfty(double new_p_infty) { p_infty = new_p_infty; };
	void calcWaveStrengths(double u, double sound, double rho, double dx, double drho, double dp, double du, StateVector2D & c, StateVector2D & cneg);
};

class SupersonicInlet :
	public Boundary
{
private:
	StateVector2D in_state;
public:
	SupersonicInlet(StateVector2D &new_in_state) : in_state(new_in_state) {};
	void setInState(StateVector2D &new_in_state) { in_state = new_in_state; };
	StateVector2D getGhostState(std::array<int, 2> ind);
};

class SlipWall :
	public Boundary
{
public:
	StateVector2D getGhostState(std::array<int, 2> ind);
};