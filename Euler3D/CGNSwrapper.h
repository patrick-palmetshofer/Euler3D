#pragma once
#include <string>
#include "Grid.h"

class CGNSwrapper
{
private:
public:
	CGNSwrapper();
	~CGNSwrapper();

	void readSingleBlockStructured(std::string filename, Eigen::Matrix<Eigen::Array2d, -1, -1>& points);
};

