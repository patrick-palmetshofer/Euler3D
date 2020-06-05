#pragma once
#include "Reconstruct.h"
#include "GlobalTypes.h"
#include "Limiter.h"
#include "Grid.h"

#include <Eigen/Dense>

class ReconstructMUSCL :
	public Reconstruct
{
private:
	std::unique_ptr<Limiter> limiter;

	//MUSCL parameters
	double reconstruct_eps = 1;
	double reconstruct_kappa = -1;

	StateVector leftdiff, diff, rightdiff;
	StateVector left_diff, right_diff;
	StateVector slope_left, slope_right;
public:
	ReconstructMUSCL();
	~ReconstructMUSCL();

	void setKappa(double kappa) { reconstruct_kappa = kappa; };
	void setEps(double eps) { reconstruct_eps = eps; };

	void setLimiter(std::unique_ptr<Limiter>& new_limiter) { limiter = std::move(new_limiter); };
	std::pair<StateVector, StateVector> reconstructStatesMUSCL(const StateVector& c_left_left, const StateVector& c_left, const StateVector& c_right, const StateVector& c_right_right);

	std::pair<StateVector, StateVector> reconstructStates(int i, int j, int k, Eigen::Array3i &dim_arr);
};

