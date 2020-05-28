#include "ReconstructMUSCL.h"
#include <Eigen/Dense>


ReconstructMUSCL::ReconstructMUSCL()
{
}


ReconstructMUSCL::~ReconstructMUSCL()
{
}

//std::pair<StateVector2D, StateVector2D> ReconstructMUSCL::reconstructStates(const StateVector2D& c_left_left, const StateVector2D& c_left, const StateVector2D& c_right, const StateVector2D& c_right_right, Eigen::Vector3d n)
//{
//	leftdiff = (c_left - c_left_left) * 2.0 / (Sx_left_left + Sx_left);
//	diff = (c_right - c_left) * 2.0 / (Sx_left + Sx_right);
//	rightdiff = (c_right_right - c_right) * 2.0 / (Sx_right_right + Sx_right);
//
//	slope_left = 0.25*reconstruct_eps*((1 + reconstruct_kappa)*diff + (1 - reconstruct_kappa)*leftdiff);
//	slope_right = 0.25*reconstruct_eps*((1 + reconstruct_kappa)*diff + (1 - reconstruct_kappa)*rightdiff);
//
//	slope_left = slope_left * limiter->limiter(leftdiff, diff);
//	slope_right = slope_right * limiter->limiter(diff, rightdiff);
//
//	reconstruct_left = c_left + slope_left * Sx_left;
//	reconstruct_right = c_right - slope_right * Sx_right;
//
//	return std::make_pair<StateVector2D, StateVector2D>(reconstruct_left, reconstruct_right);
//}

std::pair<StateVector2D, StateVector2D> ReconstructMUSCL::reconstructStatesMUSCL(const StateVector2D& c_left_left, const StateVector2D& c_left, const StateVector2D& c_right, const StateVector2D& c_right_right)
{
	leftdiff = (c_left - c_left_left);
	diff = (c_right - c_left);
	rightdiff = (c_right_right - c_right);

	slope_left = 0.25*reconstruct_eps*((1 + reconstruct_kappa)*diff + (1 - reconstruct_kappa)*leftdiff);
	slope_right = 0.25*reconstruct_eps*((1 + reconstruct_kappa)*diff + (1 - reconstruct_kappa)*rightdiff);

	slope_left = slope_left * limiter->limiter(leftdiff, diff);
	slope_right = slope_right * limiter->limiter(diff, rightdiff);

	reconstruct_left = c_left + slope_left;
	reconstruct_right = c_right - slope_right;

	return std::make_pair(reconstruct_left, reconstruct_right);
}

//Higher-order flux calculation. Reconstructs states and then calls first order flux
std::pair<StateVector2D, StateVector2D> ReconstructMUSCL::reconstructStates(int i, int j, int dim)
{
	int ineg = i, ipos = i, inegneg = i;
	int jneg = j, jpos = j, jnegneg = j;
	bool fallback = false;
	switch (dim)
	{
	case 0:
		ineg--;
		ipos++;
		inegneg -= 2;
		break;
	case 1:
		jneg--;
		jpos++;
		jnegneg -= 2;
		break;
	default:
		throw;
	}
	//Fall back if near boundary
	if (inegneg < 0 || jnegneg < 0 || jpos >= grid->getnetaCells() || ipos >= grid->getnxiCells())
		return Reconstruct::reconstructStates(i, j, dim);
	return reconstructStatesMUSCL((*conservative)(inegneg,jnegneg), (*conservative)(ineg,jneg), (*conservative)(i,j), (*conservative)(ipos,jpos));
}