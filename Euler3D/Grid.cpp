#include "Grid.h"

#include <fstream>
#include <stdexcept>
#include <iostream>

#include "CGNSwrapper.h"

#include <CGNS/src/cgnslib.h>

Grid::Grid()
{
}


Grid::~Grid()
{
}

void Grid::readCGNS(std::string filename)
{
	CGNSwrapper cgns;
	cgns.readSingleBlockStructured(filename,points);
	ncomp_cells[0] = (int) points.size()-1;
	ncomp_cells[1] = (int) points(0).size()-1;
	ncomp_cells[2] = (int) points(0)(0).size()-1;
	allocate();
}

//Reads a GridPro grid file. Structured grid.
void Grid::readGridPro(std::string filename)
{
	std::ifstream stream;

	//Number of points in both computational directions
	IndArray n_points;
	//std::vector<std::vector<Eigen::Array2D>> points;

	try
	{
		stream.open(filename, std::ifstream::in);
		if (stream.fail())
			throw;
		stream >> n_points[1];
		stream >> n_points[0];
		stream >> n_points[2];

		ncomp_cells = n_points - 1;

		std::vector<std::vector<std::vector<std::array<double, 3>>>> testpoints;

		Euler::resize(testpoints, n_points);

		Euler::resize(points,n_points);
		for (int j = n_points[1] - 1; j >= 0; j--)
		{
			for (int i = n_points[0] - 1; i >= 0; i--)
			{
				for (int k = n_points[2] - 1; k >= 0; k--)
				{
					stream >> points(i)(j)(k)[0];
					stream >> points(i)(j)(k)[1];
					stream >> points(i)(j)(k)[2];

					testpoints[i][j][k][0] = points(i)(j)(k)[0];
					testpoints[i][j][k][1] = points(i)(j)(k)[1];
					testpoints[i][j][k][2] = points(i)(j)(k)[2];
				}
			}
		}
		//for (int j = 0; j < n_points[1]; j++)
		//{
		//	for (int i = 0; i < n_points[0]; i++)
		//	{
		//		for (int k = 0; k < n_points[2]; k++)
		//		{
		//			stream >> points(i)(j)(k)[0];
		//			stream >> points(i)(j)(k)[1];
		//			stream >> points(i)(j)(k)[2];
		//		}
		//	}
		//}
		stream.close();
	}
	catch (std::ifstream::failure e) {
		std::cerr << "Exception reading file\n";
	}
	allocate();
}

void Grid::allocate()
{
	//Allocate all matrices, Allocates 1 too much in each direction. Possibly better nonethless due to data locality
	ncomp_phys.resize(dims);
	for (int dim = 0; dim < dims; ++dim)
	{
		//Fill matrices
		IndArray max_arr = ncomp_cells;
		++max_arr[dim];

		Euler::resize(ncomp_phys(dim), max_arr);
		Euler::resize(face_areas(dim), max_arr);

		std::array<DirVector, 2> face_vec;

		for (int i = 0; i < max_arr[0]; i++)
		{
			for (int j = 0; j < max_arr[1]; j++)
			{
				for (int k = 0; k < max_arr[2]; k++)
				{
					IndArray ind = { i,j,k };

					//Calculate Face Parallel Vectors
					int n = 0;
					for (int m = 0; m < dims; ++m)
					{
						if (m != dim) //Only use elements which are not the computational element
						{
							IndArray indplus = ind;
							indplus[m] = ind[m] + 1;
							auto plusp = points(indplus[0])(indplus[1])(indplus[2]);
							auto p = points(i)(j)(k);
							face_vec[n] = plusp - p;
							n++;
						}
					}

					auto Scomp_phys = face_vec[0].cross(face_vec[1]);
					//if (dim == 1)
					//	Scomp_phys = face_vec[1].cross(face_vec[0]);
					double Scomp = Scomp_phys.norm();

					if (Scomp == 0)
						throw;

					ncomp_phys(dim)(i)(j)(k) = Scomp_phys / Scomp;

					face_areas(dim)(i)(j)(k) = Scomp;
				}
			}
		}
	}

	Eigen::Matrix3d mat;

	Euler::resize(volumes, ncomp_cells);
	for (int i = 0; i < ncomp_cells[0]; ++i)
	{
		for (int j = 0; j < ncomp_cells[1]; ++j)
		{
			for (int k = 0; k < ncomp_cells[2]; ++k)
			{
				mat.col(0) = (points(i + 1)(j)(k) - points(i)(j)(k));
				mat.col(1) = (points(i)(j + 1)(k) - points(i)(j)(k));
				mat.col(2) = (points(i)(j)(k + 1) - points(i)(j)(k));
				volumes[i][j][k] = std::abs(mat.determinant());
			}
		}
	}

	DirVector max, min;
	max.fill(0);
	min.fill(1e9);
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = 0; j < points(0).size(); ++j)
		{
			for (int k = 0; k < points(0)(0).size(); ++k)
			{
				for (int dim = 0; dim < DirVector().size(); ++dim)
				{
					max[dim] = std::max(max[dim], points(i)(j)(k)[dim]);
					min[dim] = std::min(min[dim], points(i)(j)(k)[dim]);
				}
			}
		}
	}
	maxDistance = max - min;
}