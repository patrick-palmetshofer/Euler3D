#include "CGNSwrapper.h"

#include <CGNS/src/cgnslib.h>

CGNSwrapper::CGNSwrapper()
{
}


CGNSwrapper::~CGNSwrapper()
{
}

void CGNSwrapper::readSingleBlockStructured(std::string filename, GridTensor<DirVector> &points)
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

	IndArray irmaxind;
	for (int i = 0; i < irmaxind.size(); ++i)
		irmaxind[i] = (int)irmax[irmaxind.size() - 1 - i];
	Euler::resize(points, irmaxind);

	std::vector<std::vector<std::vector<std::array<double, 3>>>> testpoints;
	Euler::resize(testpoints, irmaxind);

	for (int i = 0; i < irmaxind[0]; ++i)
	{
		for (int j = 0; j < irmaxind[1]; ++j)
		{
			for (int k = 0; k < irmaxind[2]; ++k)
			{
				points(i)(j)(k)[0] = x[i*irmax[0] * irmax[1] + j * irmax[0] + k];
				points(i)(j)(k)[1] = y[i*irmax[0] * irmax[1] + j * irmax[0] + k];
				points(i)(j)(k)[2] = z[i*irmax[0] * irmax[1] + j * irmax[0] + k];

				testpoints[i][j][k][0] = x[i*irmax[0] * irmax[1] + j * irmax[0] + k];
				testpoints[i][j][k][1] = y[i*irmax[0] * irmax[1] + j * irmax[0] + k];
				testpoints[i][j][k][2] = z[i*irmax[0] * irmax[1] + j * irmax[0] + k];
			}
		}
	}

	delete[] x;
	delete[] y;
	delete[] z;
}

void CGNSwrapper::modifySingleBlockStructured(std::string filename, StateTensor &solution)
{
	int file_ind, base_ind, zone_ind, flow_ind, cell_dim, phys_dim;
	cgsize_t i_size[3][3];
	cgsize_t irmin[3];
	int irmax[3];

	char buf[256];
	if (cg_open(filename.c_str(), CG_MODE_MODIFY, &file_ind))
		throw;
	//Use only one base and zone
	base_ind = 1;
	zone_ind = 1;
	flow_ind = 1;

	cg_sol_write(file_ind, base_ind, zone_ind, "Solution1", CG_CellCenter, &flow_ind);

	//set min indices
	for (int i = 0; i < 3; ++i)
		irmin[i] = 1;

	//double x[100][100];
	//double y[100][100];
	irmax[0] = solution.size();
	irmax[1] = solution(0).size();
	irmax[2] = solution(0)(0).size();

	double* rho_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* rhou_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* rhov_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* rhow_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* et_ptr = new double[irmax[0] * irmax[1] * irmax[2]];


	for (int i = 0; i < irmax[0]; ++i)
	{
		for (int j = 0; j < irmax[1]; ++j)
		{
			for (int k = 0; k < irmax[2]; ++k)
			{
				rho_ptr[i*irmax[0] * irmax[1] + j * irmax[0] + k] = solution(i)(j)(k)[0];
				rhou_ptr[i*irmax[0] * irmax[1] + j * irmax[0] + k] = solution(i)(j)(k)[1];
				rhov_ptr[i*irmax[0] * irmax[1] + j * irmax[0] + k] = solution(i)(j)(k)[2];
				rhow_ptr[i*irmax[0] * irmax[1] + j * irmax[0] + k] = solution(i)(j)(k)[3];
				et_ptr[i*irmax[0] * irmax[1] + j * irmax[0] + k] = solution(i)(j)(k)[4];
			}
		}
	}

	int field_inds[5] = { 0,1,2,3,4 };

	//cg_section_read(file_ind,base_ind,zone_ind,"CoordinateX",)
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "Density", rho_ptr, field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "MomentumX", rhou_ptr, field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "MomentumY", rhov_ptr, field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "MomentumZ", rhow_ptr, field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "EnergyStagnationDensity", et_ptr, field_inds))
		throw;

	cg_close(file_ind);



	delete[] rho_ptr;
	delete[] rhou_ptr;
	delete[] rhov_ptr;
	delete[] rhow_ptr;
	delete[] et_ptr;

}

void CGNSwrapper::writeSingleBlockStructured(std::string filename, const GridTensor<DirVector> &points, StateTensor &solution)
{
	int file_ind, base_ind, zone_ind, coord_ind, flow_ind;
	cgsize_t i_size[3][3];
	cgsize_t irmin[3], irmaxp[3];
	int irmax[3];

	char buf[256];
	if (cg_open(filename.c_str(), CG_MODE_WRITE, &file_ind))
		throw;

	cg_base_write(file_ind, "Base", 3, 3, &base_ind);

	//set min indices
	for (int i = 0; i < 3; ++i)
		irmin[i] = 1;


	//double x[100][100];
	//double y[100][100];
	irmax[0] = solution.size();
	irmax[1] = solution(0).size();
	irmax[2] = solution(0)(0).size();


	for (int i = 0; i < 3; ++i)
		irmaxp[i] = irmax[i] + 1;

	for (int i = 0; i < 3; ++i)
	{
		i_size[0][i] = irmaxp[i];
	}
	for (int i = 0; i < 3; ++i)
	{
		i_size[1][i] = irmax[i];
	}
	for (int i = 0; i < 3; ++i)
	{
		i_size[2][i] = 0;
	}

	if (cg_zone_write(file_ind, base_ind, "Zone 1", *i_size, CG_Structured, &zone_ind))
		throw;

	double* x = new double[irmaxp[0] * irmaxp[1] * irmaxp[2]];
	double* y = new double[irmaxp[0] * irmaxp[1] * irmaxp[2]];
	double* z = new double[irmaxp[0] * irmaxp[1] * irmaxp[2]];

	for (int i = 0; i < irmaxp[0]; ++i)
	{
		for (int j = 0; j < irmaxp[1]; ++j)
		{
			for (int k = 0; k < irmaxp[2]; ++k)
			{
				int ind_unwrapped = k * irmaxp[0] * irmaxp[1] + j * irmaxp[0] + i;
				x[ind_unwrapped] = points(i)(j)(k)(0);
				y[ind_unwrapped] = points(i)(j)(k)(1);
				z[ind_unwrapped] = points(i)(j)(k)(2);
			}
		}
	}

	if (cg_coord_write(file_ind, base_ind, zone_ind, CG_RealDouble, "CoordinateX", x, &coord_ind))
		throw;
	if (cg_coord_write(file_ind, base_ind, zone_ind, CG_RealDouble, "CoordinateY", y, &coord_ind))
		throw;
	if (cg_coord_write(file_ind, base_ind, zone_ind, CG_RealDouble, "CoordinateZ", z, &coord_ind))
		throw;

	double* rho_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* rhou_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* rhov_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* rhow_ptr = new double[irmax[0] * irmax[1] * irmax[2]];
	double* et_ptr = new double[irmax[0] * irmax[1] * irmax[2]];


	for (int i = 0; i < irmax[0]; ++i)
	{
		for (int j = 0; j < irmax[1]; ++j)
		{
			for (int k = 0; k < irmax[2]; ++k)
			{
				int ind_unwrapped = k * irmax[0] * irmax[1] + j * irmax[0] + i;
				rho_ptr[ind_unwrapped] = solution(i)(j)(k)[0];
				rhou_ptr[ind_unwrapped] = solution(i)(j)(k)[1];
				rhov_ptr[ind_unwrapped] = solution(i)(j)(k)[2];
				rhow_ptr[ind_unwrapped] = solution(i)(j)(k)[3];
				et_ptr[ind_unwrapped] = solution(i)(j)(k)[4];
			}
		}
	}

	int field_inds = 0;

	cg_sol_write(file_ind, base_ind, zone_ind, "Solution1", CG_CellCenter, &flow_ind);

	//cg_section_read(file_ind,base_ind,zone_ind,"CoordinateX",)
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "Density", rho_ptr, &field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "MomentumX", rhou_ptr, &field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "MomentumY", rhov_ptr, &field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "MomentumZ", rhow_ptr, &field_inds))
		throw;
	if (cg_field_write(file_ind, base_ind, zone_ind, flow_ind, CG_RealDouble, "EnergyStagnationDensity", et_ptr, &field_inds))
		throw;

	cg_close(file_ind);

	delete[] x;
	delete[] y;
	delete[] z;
	delete[] rho_ptr;
	delete[] rhou_ptr;
	delete[] rhov_ptr;
	delete[] rhow_ptr;
	delete[] et_ptr;
}