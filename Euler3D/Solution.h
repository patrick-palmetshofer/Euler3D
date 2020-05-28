#pragma once
#include "GlobalTypes.h"
#include "Fluid.h"
#include <string>
class Solution
{
private:
	std::string filename;
	Fluid * fluid;
	//Total residual. Can be calculated through L2 or Linfty. See corresponding methods
	StateVector2D residuals_L2;
	StateVector2D residuals_Linfty;
	StateMatrix2D p;

	std::vector<StateVector2D> v_residuals_L2;
	std::vector<StateVector2D> v_residuals_Linfty;
public:
	Solution();
	~Solution();

	void setFilename(std::string& f) { filename = f; };
	void setFluid(Fluid * n_fluid) { fluid = n_fluid; };

	//Calculate residuals, total, not per equation in terms of old conservatives and new conservatives
	double calcResidualL2(StateMatrix2D & o, StateMatrix2D & n);
	double calcResidualLinfty(StateMatrix2D & o, StateMatrix2D & n);

	//Calculate residuals per value
	StateVector2D calcResidualsL2(StateMatrix2D & o, StateMatrix2D & n);
	StateVector2D calcResidualsLinfty(StateMatrix2D & o, StateMatrix2D & n);

	StateVector2D getResidualsL2();
	StateVector2D getResidualsLinfty();

	//Write solution to file. Same format as input grid
	void writeSolution(std::string filename, StateMatrix2D &c);
	void writeSolution(StateMatrix2D &c);
	void writeSolution(StateMatrix2D &c, int i);
	void writeResidual(std::vector<StateVector2D>& residuals, std::string filename);
};

