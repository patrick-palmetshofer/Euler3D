#pragma once
#include "Flux.h"
#include "Fluid.h"
class FluxRoe :
	public Flux
{
private:
	Eigen::Matrix4d roevectors;
	StateVector2D roefactors;
	StateVector2D eigenvals;
	StateVector2D dissip;
	StateVector2D prim_left;
	StateVector2D prim_right;
public:
	FluxRoe(int new_dim);
	~FluxRoe();

	//Options for numerical dissipation calculation: Roe flux (no entropy correction)
	StateVector2D calcDissip(std::pair<StateVector2D, StateVector2D> leftrightstates, double nx, double ny);
};

