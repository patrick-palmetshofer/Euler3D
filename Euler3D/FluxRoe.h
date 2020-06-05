#pragma once
#include "Flux.h"
#include "Fluid.h"
class FluxRoe :
	public Flux
{
private:
	Eigen::Matrix<double,5,5> roevectors;
	StateVector roefactors;
	StateVector eigenvals;
	StateVector dissip;
	StateVector prim_left;
	StateVector prim_right;
public:
	FluxRoe(int new_dim);
	~FluxRoe();

	//Options for numerical dissipation calculation: Roe flux (no entropy correction)
	StateVector calcDissip(std::pair<StateVector, StateVector> leftrightstates, DirVector &n);
};

