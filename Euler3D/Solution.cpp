#include "Solution.h"
#include <fstream>
#include <iostream>

#include "CGNSwrapper.h"

#ifdef _WIN32
#include <io.h> 
#define access    _access_s
#else
#include <unistd.h>
#endif

bool FileExists(const std::string &filename)
{
	return access(filename.c_str(), 0) == 0;
}

Solution::Solution()
{
}


Solution::~Solution()
{
}

//Residuals calculation. Either returning the maximum residual or residual in each variable.
StateVector Solution::calcResidualsL2(StateTensor &o, StateTensor &n)
{
	StateVector sumres;
	for (int i = 0; i < n.size(); i++)
	{
		for (int j = 0; j < n(0).size(); j++)
		{
			for (int k = 0; k < n(0)(0).size(); k++)
			{
				auto temp = (n(i)(j)(k) - o(i)(j)(k)) / ref_state.array();
				sumres += temp.matrix().squaredNorm();
			}
		}
	}

	sumres = (sumres / (n.size()*n(0).size()*n(0)(0).size())).sqrt();
	v_residuals_L2.push_back(sumres);
	return sumres;
}

StateVector Solution::calcResidualsLinfty(StateTensor &o, StateTensor &n)
{
	StateVector res;
	for (int i = 0; i < n.size(); i++)
	{
		for (int j = 0; j < n(0).size(); j++)
		{
			for (int k = 0; k < n(0)(0).size(); k++)
			{
				auto temp = (n(i)(j)(k) - o(i)(j)(k)) / ref_state.array();
				for (int m = 0; m < n(i)(j)(k).size(); ++m)
				{
					if (temp[m] > res[m])
						res[m] = temp[m];
				}
			}
		}
	}
	res = (res / (n.size()*n(0).size()*n(0)(0).size())).sqrt();
	v_residuals_L2.push_back(res);
	return res;
}


StateVector Solution::getResidualsL2()
{
	return residuals_L2;
}

StateVector Solution::getResidualsLinfty()
{
	return residuals_Linfty;
}

void Solution::writeRaw(std::string filename, StateTensor &c)
{
	//Write a set of primitive variables, with Temperature in place of e_t
	if (p.size() != c.size())
		p = StateTensor(c);
	for (int i = 0; i < c.size(); i++)
	{
		for (int j = 0; j < c(0).size(); j++)
		{
			for (int k = 0; k < c(0)(0).size(); k++)
			{
				p(i)(j)(k) = fluid->cons2prim(c(i)(j)(k));
				p(i)(j)(k).tail(1) = fluid->getTemperature(c(i)(j)(k));
			}
		}
	}

	//Write to file
	std::ofstream stream;
	try
	{
		stream.open(filename + ".sol", std::ofstream::out);
		stream << p.size() << ",\t" << p(0).size() << ",\t" << p(0)(0).size() << "\n";
		stream << "rho" << ",\t" << "u" << ",\t" << "v" << ",\t" << "w" << ",\t" << "T" << "\n";
		//Reserve add here
		for (int i = 0; i < p.size(); i++) 
		{
			for (int j = 0; j < p(0).size(); j++) 
			{
				for (int k = 0; k < p(0)(0).size(); k++) 
				{
					for (int dim = 0; dim < p(0)(0)(0).size() - 1; ++dim)
						stream << p(i)(j)(k)[dim] << ",\t";
					stream << p(i)(j)(k)[p(0)(0)(0).size() - 1] << "\n";
				}
			}
		}
		stream.flush();
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}

void writeCGNS(std::string filename, Grid * grid, StateTensor &c)
{
	CGNSwrapper cgns;
	cgns.writeSingleBlockStructured(filename, grid->getPoints(), c);
}

void writeCGNS(std::string filename, StateTensor &c)
{
	CGNSwrapper cgns;
	cgns.modifySingleBlockStructured(filename, c);		
}

void Solution::writeSolution(std::string filename, Grid * grid, StateTensor &c)
{
	writeCGNS(filename, grid, c);
	writeResidual(v_residuals_L2, filename + "_L2.res");
	writeResidual(v_residuals_Linfty, filename + "_Linfty.res");
}

void Solution::writeSolution(std::string filename, StateTensor &c)
{
	writeCGNS(filename, c);
	writeResidual(v_residuals_L2, filename + "_L2.res");
	writeResidual(v_residuals_Linfty, filename + "_Linfty.res");
}

void Solution::writeSolution(StateTensor & c)
{
	writeSolution(filename, c);
}

void Solution::writeSolution(StateTensor & c, int i)
{
	writeSolution(filename+std::to_string(i), c);
}

void Solution::writeResidual(std::vector<StateVector> &residuals, std::string filename)
{
	std::ofstream stream;
	try
	{
		stream.open(filename, std::ofstream::out);
		stream << "rho" << ",\t" << "rhou" << ",\t" << "rhov" << ",\t" << "rhoet" << "\n";
		//Reserve add here
		for (int i = 0; i < residuals.size(); i++)
		{
			for (int k = 0; k < 3; k++)
			{
				stream << residuals[i][k] << ",\t";
			}
			stream << residuals[i][3] << "\n";
		}
		stream.flush();
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}