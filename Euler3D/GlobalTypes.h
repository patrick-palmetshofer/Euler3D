#pragma once
#include <memory>
#include <vector>
#include <valarray>
#include <array>

#include <iostream>

#include <Eigen/Dense>

//Generalized State Vector specified to 2 dimensions
using StateVector = Eigen::Array<double,5,1>;
using DirVector = Eigen::Vector3d;
using IndArray = Eigen::Array3i;

template <typename T>
using GridTensor = Eigen::Array<Eigen::Array<Eigen::Array<T, -1, 1>, -1, 1>, -1, 1>;

//Matrices of StateVectors or Values defined on the grid
using ValueTensor = GridTensor<double>;
using StateTensor = GridTensor<StateVector>;
//using Matrix2D = std::vector<std::vector<std::vector<stateVector2D>>>;

//template <uint ndim>
//using Jacobian = std::array<std::array<double,ndim+2>,ndim+2>;

//using Jacobian2D = Jacobian<2>;

namespace Euler
{
	template<typename T>
	void resize(T &grid, IndArray max_inds)
	{
		grid.resize(max_inds[0]);
		for (int i = 0; i < max_inds[0]; i++)
		{
			grid[i].resize(max_inds[1]);
			for (int j = 0; j < max_inds[1]; j++)
			{
				grid[i][j].resize(max_inds[2]);
			}
		}
	}

	template<typename T>
	void fill(GridTensor<T> &grid, T value)
	{
		for (int i = 0; i < grid.size(); i++)
		{
			for (int j = 0; j < grid(0).size(); j++)
			{
				for (int k = 0; k < grid(0)(0).size(); k++)
				{
					grid(i)(j)(k) = value;
				}
			}
		}
	}

	//Utilities to swap coordinates in Statevectors
	template<typename T>
	T swap(T &data, int ind1, int ind2)
	{
		if (ind1 == ind2)
			return data;
		T res = data;
		res[ind1] = data[ind2];
		res[ind2] = data[ind1];
		return res;
	}

	template<typename T>
	T swap(T &data, int i)
	{
		return Euler::swap(data, 1, i);
	}


	inline bool checkNaN(StateVector *m)
	{
		for (int k = 0; k < m->size(); k++)
		{
			if (!std::isfinite((*m)[k]))
			{
				return true;
			}
		}
		return false;
	}

	//Checks for errors in matrices, used for debugging
	inline bool checkNaN(StateTensor *m)
	{
		size_t xi_size = m->size();
		size_t eta_size = (*m)(0).size();
		size_t zeta_size = (*m)(0)(0).size();
		for (int i = 0; i < xi_size; i++)
		{
			for (int j = 0; j < eta_size; j++)
			{
				for (int k = 0; k < zeta_size; k++)
				{
					if (checkNaN(&(*m)(i)(j)(k)))
						return true;
				}
			}
		}
		return false;
	}
};


//namespace std {
//	template<class T> struct _Unique_if {
//		typedef unique_ptr<T> _Single_object;
//	};
//
//	template<class T> struct _Unique_if<T[]> {
//		typedef unique_ptr<T[]> _Unknown_bound;
//	};
//
//	template<class T, size_t N> struct _Unique_if<T[N]> {
//		typedef void _Known_bound;
//	};
//
//	template<class T, class... Args>
//	typename _Unique_if<T>::_Single_object
//		make_unique(Args&&... args) {
//		return unique_ptr<T>(new T(std::forward<Args>(args)...));
//	}
//
//	template<class T>
//	typename _Unique_if<T>::_Unknown_bound
//		make_unique(size_t n) {
//		typedef typename remove_extent<T>::type U;
//		return unique_ptr<T>(new U[n]());
//	}
//
//	template<class T, class... Args>
//	typename _Unique_if<T>::_Known_bound
//		make_unique(Args&&...) = delete;
//}