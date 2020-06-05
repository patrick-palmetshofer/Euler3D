#include "ReconstructMUSCL.h"
#include <Eigen/Dense>


ReconstructMUSCL::ReconstructMUSCL()
{
}


ReconstructMUSCL::~ReconstructMUSCL()
{
}

//std::pair<StateVector, StateVector> ReconstructMUSCL::reconstructStates(const StateVector& c_left_left, const StateVector& c_left, const StateVector& c_right, const StateVector& c_right_right, Eigen::Vector3d n)
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
//	return std::make_pair<StateVector, StateVector>(reconstruct_left, reconstruct_right);
//}

std::pair<StateVector, StateVector> ReconstructMUSCL::reconstructStatesMUSCL(const StateVector& c_left_left, const StateVector& c_left, const StateVector& c_right, const StateVector& c_right_right)
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
std::pair<StateVector, StateVector> ReconstructMUSCL::reconstructStates(int i, int j, int k, IndArray &dim_arr)
{
	IndArray ijk = { i,j,k };
	IndArray ijknegneg = ijk - 2 * dim_arr;
	IndArray ijkneg = ijk - dim_arr;
	IndArray ijkpos = ijk + dim_arr;

	//Fall back if near boundary
	for (int dim = 0; dim < dim_arr.size(); ++dim)
	{
		if (ijknegneg[dim] < 0 || ijkpos[dim] >= grid->getnComponentCells(dim))
			return Reconstruct::reconstructStates(i, j, k, dim_arr);
	}
	return reconstructStatesMUSCL((*conservative)(ijknegneg[0])(ijknegneg[1])(ijknegneg[2]), (*conservative)(ijkneg[0])(ijkneg[1])(ijkneg[2]), (*conservative)(i)(j)(k), (*conservative)(ijkpos[0])(ijkpos[1])(ijkpos[2]));
}