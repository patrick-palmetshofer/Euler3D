#pragma once
#include "GlobalTypes.h"

class Grid
{
private:
	//Numbers of cells in both computational directions, excluding ghost cells
	int nxi_cells;
	int neta_cells;

	//2D Vector of points. 
	Eigen::Matrix<Eigen::Array2d,-1,-1> points;

	Eigen::Array2d maxDistance;

	//Face normal vector components for each computational direction (2 physical directions each)
	ValueMatrix2D nxi_xs;
	ValueMatrix2D nxi_ys;
	ValueMatrix2D neta_xs;
	ValueMatrix2D neta_ys;

	//Face areas in both computational directions
	ValueMatrix2D Sxis;
	ValueMatrix2D Setas;

	//Volumes of cells in domain
	ValueMatrix2D volumes;
public:
	Grid();
	~Grid();

	Eigen::Array2d getPoint(int i, int j) { return points(i,j);	};

	double getMaxDistance(int dim) { return maxDistance[dim]; };

	int getnxiCells() { return nxi_cells; };
	int getnetaCells() { return neta_cells; };

	int getnComponentCells(int dim) { return dim ? neta_cells : nxi_cells; };

	double getnXiXs(int i, int j) { return nxi_xs(i,j); };
	double getnEtaXs(int i, int j) { return neta_xs(i,j); };
	double getnXiYs(int i, int j) { return nxi_ys(i,j); };
	double getnEtaYs(int i, int j) { return neta_ys(i,j); };

	double getSxi(int i, int j) { return Sxis(i,j); };
	double getSeta(int i, int j) { return Setas(i,j); };

	double getVolume(int i, int j) { return volumes(i,j); };

	double getnComponent(int i, int j, int compdim, int physdim);

	void readCGNS(std::string filename);

	//Reads GridPro grid file (Warning, no treatment of Gridpro Boundary Conditions)
	void readGridPro(std::string filename);
	void allocate();
};

