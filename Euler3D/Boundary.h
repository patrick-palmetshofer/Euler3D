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
	StateTensor * conservative;
	Flux * flux;

	int dim;
	int dir;
	IndArray dir_arr;
	std::vector<IndArray> cell_inds;

	StateVector ghost_state;

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

	void setFront();

	void setBack();

	void setDimDir(int dim, int dir);

	void init();

	void setGrid(Grid * new_grid) { grid = new_grid; };
	void setFluid(Fluid * newfluid) { fluid = newfluid; };
	void setFlux(Flux * newflux) { flux = newflux; };
	void setConservative(StateTensor * cons) { conservative = cons; };

	int getDim() { return dim; };
	int getDir() { return dir; };

	virtual StateVector getGhostState(IndArray &ind) = 0;

	void apply();
};

class BoundaryLODI :
	public Boundary
{
protected:
	double L1 = 0, L2 = 0, L3 = 0, L4 = 0, L5 = 0;
public:
	virtual void calcWaveStrengths(double u, double sound, double rho, double dx, double drho, double dp, double du, StateVector & c, StateVector & cneg) = 0;
	StateVector getGhostState(IndArray &ind);
};

class SupersonicOutlet :
	public Boundary
{
public:
	//StateVector getGhostState(std::array<int, 2> ind);
	StateVector getGhostState(IndArray& ind);
};

class Outlet :
	public BoundaryLODI
{
private:
	double p_infty;
public:
	Outlet(double new_pinfty) : p_infty(new_pinfty) {};
	void setPinfty(double new_p_infty) { p_infty = new_p_infty; };
	void calcWaveStrengths(double u, double sound, double rho, double dx, double drho, double dp, double du, StateVector & c, StateVector & cneg);
};

class SupersonicInlet :
	public Boundary
{
private:
	StateVector in_state;
public:
	SupersonicInlet(StateVector &new_in_state) : in_state(new_in_state) {};
	void setInState(StateVector &new_in_state) { in_state = new_in_state; };
	StateVector getGhostState(IndArray& ind);
};

class SlipWall :
	public Boundary
{
public:
	StateVector getGhostState(IndArray& ind);
};