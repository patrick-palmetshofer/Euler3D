#pragma once
#include <string>
#include "Grid.h"

class CGNSwrapper
{
private:
public:
	CGNSwrapper();
	~CGNSwrapper();

	void readSingleBlockStructured(std::string filename, GridTensor<DirVector>& points);
	void modifySingleBlockStructured(std::string filename, StateTensor & solution);
	void modifySingleBlockStructured(std::string filename, GridTensor<DirVector>& points, StateTensor & solution);
	void writeSingleBlockStructured(std::string filename, const GridTensor<DirVector>& points, StateTensor & solution);
};

