#pragma once
#include "Grid.h"
#include "GlobalTypes.h"

class Reconstruct
{
protected:
	StateVector2D reconstruct_left, reconstruct_right;
	Grid * grid;

	StateMatrix2D * conservative;
public:
	Reconstruct();
	virtual ~Reconstruct();

	void setConservative(StateMatrix2D * cons) { conservative = cons; };

	void setGrid(Grid * new_grid) { grid = new_grid; };

	std::pair<StateVector2D, StateVector2D> reconstructStates(const StateVector2D & c_left, const StateVector2D & c_right);

	virtual std::pair<StateVector2D, StateVector2D> reconstructStates(int i, int j, int dim);
	std::pair<StateVector2D, StateVector2D> reconstructStatesXi(int i, int j);
	std::pair<StateVector2D, StateVector2D> reconstructStatesEta(int i, int j);
};


class ReconstructFirstOrder :
	public Reconstruct
{
};