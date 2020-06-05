#pragma once
#include "GlobalTypes.h"
#include "Fluid.h"
#include "Grid.h"

#include <string>
class Solution
{
private:
	std::string filename;
	Fluid * fluid;
	//Total residual. Can be calculated through L2 or Linfty. See corresponding methods
	StateVector residuals_L2;
	StateVector residuals_Linfty;

	StateVector ref_state;

	StateTensor p;

	std::vector<StateVector> v_residuals_L2;
	std::vector<StateVector> v_residuals_Linfty;
public:
	Solution();
	~Solution();

	void setFilename(std::string& f) { filename = f; };
	void setFluid(Fluid * n_fluid) { fluid = n_fluid; };

	//Calculate residuals, total, not per equation in terms of old conservatives and new conservatives
	double calcResidualL2(StateTensor & o, StateTensor & n);
	double calcResidualLinfty(StateTensor & o, StateTensor & n);

	//Calculate residuals per value
	StateVector calcResidualsL2(StateTensor & o, StateTensor & n);
	StateVector calcResidualsLinfty(StateTensor & o, StateTensor & n);

	StateVector getResidualsL2();
	StateVector getResidualsLinfty();

	void writeRaw(std::string filename, StateTensor & c);

	void writeSolution(std::string filename, Grid * grid, StateTensor & c);

	//Write solution to file. Same format as input grid
	void writeSolution(std::string filename, StateTensor &c);
	void writeSolution(StateTensor &c);
	void writeSolution(StateTensor &c, int i);
	void writeResidual(std::vector<StateVector>& residuals, std::string filename);
};

