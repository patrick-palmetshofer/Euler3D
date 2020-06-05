#pragma once
#include "GlobalTypes.h"

class Grid
{
private:
	const int dims = 3;
	//Numbers of cells in both computational directions, excluding ghost cells
	IndArray ncomp_cells;

	//2D Vector of points. 
	GridTensor<DirVector> points;

	DirVector maxDistance;

	//Face normal vector components for each computational direction (2 physical directions each)
	Eigen::Array<GridTensor<Eigen::Vector3d>, 3, 1> ncomp_phys;

	//Face areas in both computational directions
	Eigen::Array<ValueTensor,3,1> face_areas;

	//Volumes of cells in domain
	ValueTensor volumes;
public:
	Grid();
	~Grid();

	DirVector getPoint(int i, int j, int k) { return points(i)(j)(k);	};
	DirVector getPoint(IndArray &ind) { return points(ind[0])(ind[1])(ind[2]); };
	const GridTensor<DirVector>& getPoints() { return points; };

	double getMaxDistance(int dim) { return maxDistance[dim]; };

	Eigen::Array3i getnComponentCells() { return ncomp_cells; };
	int getnComponentCells(int dim) { 
		return ncomp_cells(dim); };

	double getFaceArea(int i, int j, int k, int dim) { return face_areas(dim)(i)(j)(k); };
	double getFaceArea(IndArray ind, int dim) { return face_areas(dim)(ind[0])(ind[1])(ind[2]); };

	double getVolume(int i, int j, int k) { return volumes(i)(j)(k); };

	double getnComponent(int i, int j, int k, int comp, int phys) { return ncomp_phys(comp)(i)(j)(k)(phys); }
	Eigen::Vector3d getnVec(int i, int j, int k, int comp) { return ncomp_phys(comp)(i)(j)(k); }
	Eigen::Vector3d getnVec(IndArray ind, int comp) { return ncomp_phys(comp)(ind[0])(ind[1])(ind[2]); }


	//Reads GridPro grid file (Warning, no treatment of Gridpro Boundary Conditions)
	void readGridPro(std::string filename);
	void readCGNS(std::string filename);

	void allocate();
};

