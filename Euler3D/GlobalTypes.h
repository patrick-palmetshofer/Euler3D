#pragma once
#include <memory>
#include <vector>
#include <valarray>
#include <array>

#include <iostream>

#include <Eigen/Dense>

//Generalized State Vector specified to 2 dimensions
using StateVector2D = Eigen::Array4d;

//Matrices of StateVectors or Values defined on the grid
using ValueMatrix2D = Eigen::MatrixXd;
using StateMatrix2D = Eigen::Matrix<StateVector2D, -1, -1>;
//using Matrix2D = std::vector<std::vector<std::vector<stateVector2D>>>;

//template <uint ndim>
//using Jacobian = std::array<std::array<double,ndim+2>,ndim+2>;

//using Jacobian2D = Jacobian<2>;

namespace Euler
{
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


	inline bool checkNaN(StateVector2D *m)
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
	inline bool checkNaN(StateMatrix2D *m)
	{
		size_t xi_size = m->rows();
		size_t eta_size = (*m).cols();
		for (int i = 0; i < xi_size; i++)
		{
			for (int j = 0; j < eta_size; j++)
			{
				if (checkNaN(&(*m)(i,j)))
					return true;
			}
		}
		return false;
	}
}


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