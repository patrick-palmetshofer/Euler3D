#pragma once
#include <memory>
#include <string>
#include <map>

#include "Boundary.h"
#include "Grid.h"
#include "Solver.h"
#include "Reconstruct.h"

class Config
{
private:
	std::map<std::string, std::string> map;

	std::unique_ptr<Grid> grid;
	std::unique_ptr<Fluid> fluid;
	std::unique_ptr<Solver> solver;

	void meshRead(const std::vector<std::string>& strings);
	void fluidsRead(const std::vector<std::string>& strings);
	void solverRead(const std::vector<std::string>& strings);
	void boundaryRead(const std::vector<std::string>& strings);
	void initalRead(const std::vector<std::string>& strings);
	void solutionRead(const std::vector<std::string>& strings);

	void fillCategory(std::string & category, std::vector<std::string>& current_lines);

public:
	Config();
	~Config();

	bool cleanString(std::string & line);

	void readConfig(const std::string filename);

	std::unique_ptr<Solver> getSolver() { return std::move(solver); };
};
