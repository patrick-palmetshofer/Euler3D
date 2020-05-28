#include "Reconstruct.h"



Reconstruct::Reconstruct()
{
}


Reconstruct::~Reconstruct()
{
}

std::pair<StateVector2D, StateVector2D> Reconstruct::reconstructStates(const StateVector2D& c_left, const StateVector2D& c_right)
{
	return std::make_pair(c_left, c_right);
}

std::pair<StateVector2D, StateVector2D> Reconstruct::reconstructStates(int i, int j, int dim)
{
	int ineg = i;
	int jneg = j;
	if (dim == 0)
	{
		ineg--;
	}
	else if (dim == 1)
	{
		jneg--;
	}
	else
		throw;
	return reconstructStates((*conservative)(ineg,jneg), (*conservative)(i,j));
}

std::pair<StateVector2D, StateVector2D> Reconstruct::reconstructStatesXi(int i, int j)
{
	return reconstructStates(i, j, 0);
}

std::pair<StateVector2D, StateVector2D> Reconstruct::reconstructStatesEta(int i, int j)
{
	return reconstructStates(i, j, 1);
}
