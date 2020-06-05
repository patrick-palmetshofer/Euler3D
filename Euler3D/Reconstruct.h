#pragma once
#include "Grid.h"
#include "GlobalTypes.h"

class Reconstruct
{
protected:
	StateVector reconstruct_left, reconstruct_right;
	Grid * grid;

	StateTensor * conservative;
public:
	Reconstruct();
	virtual ~Reconstruct();

	void setConservative(StateTensor * cons) { conservative = cons; };

	void setGrid(Grid * new_grid) { grid = new_grid; };

	std::pair<StateVector, StateVector> reconstructStates(const StateVector & c_left, const StateVector & c_right);

	virtual std::pair<StateVector, StateVector> reconstructStates(int i, int j, int k, Eigen::Array3i& dim_arr);
};


class ReconstructFirstOrder :
	public Reconstruct
{
};