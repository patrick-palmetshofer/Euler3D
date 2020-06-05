#include "Reconstruct.h"



Reconstruct::Reconstruct()
{
}


Reconstruct::~Reconstruct()
{
}

std::pair<StateVector, StateVector> Reconstruct::reconstructStates(const StateVector& c_left, const StateVector& c_right)
{
	return std::make_pair(c_left, c_right);
}

std::pair<StateVector, StateVector> Reconstruct::reconstructStates(int i, int j, int k, Eigen::Array3i& dim_arr)
{
	Eigen::Array3i ijkneg = { i,j,k };
	ijkneg -= dim_arr;
	return reconstructStates((*conservative)(ijkneg[0])(ijkneg[1])(ijkneg[2]), (*conservative)(i)(j)(k));
}