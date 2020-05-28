#include "Config.h"
#include "PerfectGas.h"
#include "GlobalTypes.h"

#include <fstream>
#include <algorithm>
#include <exception>
#include <memory>
#include <cctype>

#include "TimeStepperGlobal.h"
#include "TimeStepperLocal.h"
#include "Limiter.h"
#include "ReconstructMUSCL.h"
#include "FluxRoe.h"

Config::Config()
{
	solver = std::make_unique<Solver>();
	grid = std::make_unique<Grid>();
}


Config::~Config()
{

}

bool Config::cleanString(std::string& line)
{
	std::size_t pos = line.find("//");      // position of "live" in str
	if (pos)
	{
		line = line.substr(0, pos);
	}
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
	line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
	if (line.empty())
	{
		return false;
	}
	std::transform(line.begin(), line.end(), line.begin(), [](unsigned char c) { return std::tolower(c); });
	return true;
}

void Config::meshRead(const std::vector<std::string>& strings)
{
	std::string fileext, filename;
	fileext.reserve(128);
	filename.reserve(128);
	for (auto &s : strings)
	{
		std::size_t pos = s.find("=");
		if (pos)
		{
			std::string key = s.substr(0,pos);
			std::string val = s.substr(pos+1);

			if (key == "type")
			{
				if (val != "structured" && val != "2")
					throw std::invalid_argument(val); //Only structured body-fitted meshes implemented
			}
			else if (key == "format")
			{
				fileext = val;
			}
			else if (key == "file")
			{
				filename = val;
			}
		}
	}

	if (fileext == ".grd")
		grid->readGridPro(filename + fileext);
	else if (fileext == ".cgns")
		grid->readCGNS(filename + fileext);
	else
		throw std::invalid_argument(fileext); //Only GridPro files .grd implemented
}

void Config::fluidsRead(const std::vector<std::string>& strings)
{
	for (auto &s : strings)
	{
		std::size_t pos = s.find("=");
		if (pos)
		{
			std::string key = s.substr(0, pos);
			std::string val = s.substr(pos+1);

			if (key == "type")
			{
				if (val == "perfectgas")
					fluid = std::make_unique<PerfectGas>(1.0,1.0);
			}
			else if (key == "gamma")
			{
				fluid->setGamma(std::stod(val));
			}
			else if (key == "cp")
			{
				fluid->setCp(std::stod(val));
			}
		}
	}
	fluid->setAllFromGammaCp();
}

void Config::solverRead(const std::vector<std::string>& strings)
{
	for (auto sref = strings.begin(); sref < strings.end(); ++sref)
	{
		auto &s = *sref;
		std::size_t pos = s.find("=");
		if (!pos)
		{
			continue;
		}
		std::string key = s.substr(0, pos);
		std::string val = s.substr(pos+1);

		if (key == "maxiter")
			solver->setMaxIter(std::stod(val));
		else if (key == "tend")
			solver->setMaxTime(std::stod(val));
		else if (key == "cfl")
			solver->setCFL(std::stod(val));
		else if (key == "timestepping")
		{
			std::unique_ptr<TimeStepper> stepper;
			if (val == "global")
				stepper = std::make_unique<TimeStepperGlobal>();
			else if (val == "local")
				stepper = std::make_unique<TimeStepperLocal>();
			solver->setTimeStepper(stepper);
		}
		else if (key == "implicit")
			solver->setImplicit(std::stod(val));
		else if (key == "flux")
		{
			std::unique_ptr<Flux> xi_flux, eta_flux;
			if (val == "roe:")
			{
				xi_flux = std::make_unique<FluxRoe>(0);
				eta_flux = std::make_unique<FluxRoe>(1);
			}
			solver->setXiFluxes(xi_flux);
			solver->setEtaFluxes(eta_flux);
		}
		else if (key == "reconstruct")
		{
			std::unique_ptr<Reconstruct> reconstruct;
			if (val == "firstorder" || val == "no" || val == "0" || val == "f")
				reconstruct = std::make_unique<ReconstructFirstOrder>();
			else if (val == "muscl:")
			{
				std::unique_ptr<ReconstructMUSCL> recMUSCL = std::make_unique<ReconstructMUSCL>();
				std::unique_ptr<Limiter> limiter;
				for (auto sref2 = sref; sref2 < strings.end(); ++sref2)
				{
					auto &s = *sref2;
					std::size_t pos = s.find("=");
					if (!pos)
					{
						continue;
					}
					std::string key = s.substr(0, pos);
					std::string val = s.substr(pos+1);
					if (key == "muscllimiter")
					{
						if (val == "minmod")
							limiter = std::make_unique<LimiterMinmod>();
						else if (val == "monotonecentered" || val == "mc")
							limiter = std::make_unique<LimiterMonotoneCentered>();
						else if (val == "vanalbada")
							limiter = std::make_unique<LimiterVanAlbada>();
					}
					else if (key == "musclkappa")
						recMUSCL->setKappa(std::stod(val));
					else if (key == "musclepsilon")
						recMUSCL->setEps(std::stod(val));
				}
				if (limiter == nullptr)
					limiter = std::make_unique<NoLimiter>();
				recMUSCL->setLimiter(limiter);
				reconstruct = std::move(recMUSCL);
			}
			solver->setReconstruct(reconstruct);
		}
	}
}

template<typename T>
void assignBoundary(T sref, Boundary * b)
{
	auto &s = *(sref + 1);
	std::size_t pos = s.find("=");
	std::string key = s.substr(0, pos);
	std::string val = s.substr(pos + 1);
	if (key != "cells")
		throw;

	if (val == "bottom")
		b->setBottom();
	else if (val == "top")
		b->setTop();
	else if (val == "right")
		b->setRight();
	else if (val == "left")
		b->setLeft();
	else
		throw;
}

void Config::boundaryRead(const std::vector<std::string>& strings)
{
	int nBC = 0;
	std::string s = strings.front();
	std::size_t pos = s.find("=");
	if (!pos)
	{
		throw;
	}
	std::string key = s.substr(0, pos);
	std::string val = s.substr(pos+1);

	if (key != "nbc")
		throw std::invalid_argument(key);

	nBC = std::stoi(val);
	std::vector<std::unique_ptr<Boundary>> boundaries(nBC);
	int iBC = 0;

	for (auto sref = strings.begin()+1; sref < strings.end(); ++sref)
	{
		auto &s = *sref;
		std::size_t pos = s.find("=");
		if (!pos)
		{
			continue;
		}
		std::string key = s.substr(0, pos);
		std::string val = s.substr(pos+1);

		if (val == "slipwall:")
		{
			boundaries.at(iBC) = std::make_unique<SlipWall>();
			assignBoundary(sref, boundaries.at(iBC).get());
			iBC++;
		}
		else if (val == "supersonicoutlet:")
		{
			boundaries.at(iBC) = std::make_unique<SupersonicOutlet>();
			assignBoundary(sref, boundaries.at(iBC).get());
			iBC++;
		}
		else if (val == "outlet:" || val == "subsonicoutlet:")
		{
			auto &ps = *(sref + 2);
			std::size_t ppos = ps.find("=");
			std::string pkey = ps.substr(0, ppos);
			std::string pval = ps.substr(ppos + 1);

			if (pkey != "p")
				throw;

			double p = std::stod(pval);

			boundaries.at(iBC) = std::make_unique<Outlet>(p);
			assignBoundary(sref, boundaries.at(iBC).get());
			iBC++;
		}
		else if (val == "supersonicinlet:")
		{
			double p, u, v, T;
			for (auto i = 0; i < 4; i++)
			{
				auto &states = *(sref + 2 + i);
				std::size_t statepos = states.find('=');
				std::string statekey = states.substr(0, statepos);
				std::string stateval = states.substr(statepos + 1);

				if (statekey == "p")
					p = std::stod(stateval);
				else if (statekey == "u")
					u = std::stod(stateval);
				else if (statekey == "v")
					v = std::stod(stateval);
				else if (statekey == "t")
					T = std::stod(stateval);
			}
			StateVector2D cons = fluid->user2cons(p, u, v, T);
			boundaries.at(iBC) = std::make_unique<SupersonicInlet>(cons);
			assignBoundary(sref, boundaries.at(iBC).get());
			iBC++;
		}

	}
	for (auto &b : boundaries)
	{
		if (!b)
			throw;
	}
	solver->setBoundaries(boundaries);
}

void Config::initalRead(const std::vector<std::string>& strings)
{
	double p = 1e5, u = 0, v = 0, T = 300;
	for (auto &s : strings)
	{
		std::size_t pos = s.find("=");
		if (!pos)
		{
			continue;
		}
		std::string key = s.substr(0, pos);
		std::string val = s.substr(pos+1);

		if (key == "type")
		{
			if (val != "homogenous")
				throw std::invalid_argument(val);
		}
		else if (key == "p")
			p = std::stod(val);
		else if (key == "u")
			u = std::stod(val);
		else if (key == "v")
			v = std::stod(val);
		else if (key == "T")
			T = std::stod(val);
	}
	solver->setFluid(fluid);
	solver->setConsInitial(p, u, v, T);
}

void Config::solutionRead(const std::vector<std::string>& strings)
{
	std::unique_ptr<Solution> solution = std::make_unique<Solution>();
	for (auto &s : strings)
	{
		std::size_t pos = s.find("=");
		if (!pos)
		{
			continue;
		}
		std::string key = s.substr(0, pos);
		std::string val = s.substr(pos+1);

		if (key == "filename")
			solution->setFilename(val);
		else if (key == "writeiterinterval")
			solver->setIterWriteInterval(std::stoi(val));
		else if (key == "writetimeinterval")
			solver->setTimeWriteInterval(std::stod(val));
	}
	solver->setSolution(solution);
}

void Config::fillCategory(std::string& category, std::vector<std::string>& current_lines)
{
	if (category == "mesh:")
		meshRead(current_lines);
	else if (category == "fluids:")
		fluidsRead(current_lines);
	else if (category == "solver:")
		solverRead(current_lines);
	else if (category == "boundarycondition:")
		boundaryRead(current_lines);
	else if (category == "initialcondition:")
		initalRead(current_lines);
	else if (category == "solution:")
		solutionRead(current_lines);
}

void Config::readConfig(const std::string filename)
{
	std::ifstream file;
	file.open(filename, std::ifstream::in);
	std::string line;
	line.reserve(1000);

	std::string category = "";
	std::vector<std::string> current_lines;
	current_lines.reserve(200);

	while (std::getline(file,line))
	{
		if (!cleanString(line))
			continue;

		if (line.back() == ':' && line.find("=") == std::string::npos)
		{
			fillCategory(category, current_lines);

			category = line;
			current_lines.clear();
			continue;
		}
		current_lines.push_back(line);
	}
	fillCategory(category, current_lines);
	solver->setGrid(grid);
	solver->initSizeFromGrid();

	solver->crossPopulatePointers();
}
