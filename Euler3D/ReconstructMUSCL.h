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

	StateVector2D leftdiff, diff, rightdiff;
	StateVector2D left_diff, right_diff;
	StateVector2D slope_left, slope_right;
public:
	ReconstructMUSCL();
	~ReconstructMUSCL();

	void setKappa(double kappa) { reconstruct_kappa = kappa; };
	void setEps(double eps) { reconstruct_eps = eps; };

	void setLimiter(std::unique_ptr<Limiter>& new_limiter) { limiter = std::move(new_limiter); };
	std::pair<StateVector2D, StateVector2D> reconstructStatesMUSCL(const StateVector2D& c_left_left, const StateVector2D& c_left, const StateVector2D& c_right, const StateVector2D& c_right_right);

	std::pair<StateVector2D, StateVector2D> reconstructStates(int i, int j, int dim);
};

