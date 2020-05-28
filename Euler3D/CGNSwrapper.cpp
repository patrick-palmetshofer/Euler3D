#include "CGNSwrapper.h"

#include <CGNS/src/cgnslib.h>

CGNSwrapper::CGNSwrapper()
{
}


CGNSwrapper::~CGNSwrapper()
{
}

void CGNSwrapper::readSingleBlockStructured(std::string filename, Eigen::Matrix<Eigen::Array2d, -1, -1> &points)
{
	int file_ind, base_ind, zone_ind, cell_dim, phys_dim;
	cgsize_t i_size[3][3];
	cgsize_t irmin[3], irmax[3];

	char buf[256];
	if (cg_open(filename.c_str(), CG_MODE_READ, &file_ind))
		throw;
	//Use only one base and zone
	base_ind = 1;
	zone_ind = 1;

	cg_base_read(file_ind, base_ind, buf, &cell_dim, &phys_dim);
	if (cell_dim > 3 || phys_dim > 3 || cell_dim != phys_dim)
		throw;

	cg_zone_read(file_ind, base_ind, zone_ind, buf, *i_size);
	std::string zone_name(buf);

	//Read max indices in each direction
	for (int i = 0; i < 3; ++i)
		irmax[i] = i_size[0][i];

	//Structured grid checker
	if (i_size[2][0] != 0 || i_size[2][1] != 0 || i_size[2][2] != 0)
		throw;

	//set min indices
	for (int i = 0; i < 3; ++i)
		irmin[i] = 1;

	double* x = new double[irmax[0] * irmax[1] * irmax[2]];
	double* y = new double[irmax[0] * irmax[1] * irmax[2]];
	double* z = new double[irmax[0] * irmax[1] * irmax[2]];

	//double x[100][100];
	//double y[100][100];

	//cg_section_read(file_ind,base_ind,zone_ind,"CoordinateX",)
	if (cg_coord_read(file_ind, base_ind, zone_ind, "CoordinateX", CG_RealDouble, irmin, irmax, x))
		throw;
	if (cg_coord_read(file_ind, base_ind, zone_ind, "CoordinateY", CG_RealDouble, irmin, irmax, y))
		throw;
	if (cg_coord_read(file_ind, base_ind, zone_ind, "CoordinateZ", CG_RealDouble, irmin, irmax, z))
		throw;
	
	cg_close(file_ind);

	points.resize(irmax[0], irmax[1]);

	for (int i = 0; i < irmax[0]; ++i)
	{
		for (int j = 0; j < irmax[1]; ++j)
		{
			for (int k = 0; k < irmax[2]; ++k)
			{
				points(i, j)[0] = x[j*irmax[1] + i];
				points(i, j)[1] = y[j*irmax[1] + i];
				points(i, j)[2] = y[j*irmax[1] + i];
			}
		}
	}

	delete [] x;
	delete [] y;
	delete[] z;
}
