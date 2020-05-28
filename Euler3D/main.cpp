#include "Solver.h"
#include "Config.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>

int main(int argc, char* argv[])
{

	Config config;
	std::unique_ptr<Solver> solver;

	std::string configfile = "config.cfg";
	if (argc > 1)
	{
		std::vector<std::string> argList(argv + 1, argv + argc);
		configfile = argv[0];
	}

	config.readConfig(configfile);

	solver = config.getSolver();

	solver->solve();

	return 0;
}











//const double eps = 1;
//const double kappa = 1.0/3.0;// 1.0 / 2.0;// 1.0 / 2.0;// 1.0 / 3.0;// 1.0 / 3.0;// 1.0 / 3.0;
//const bool limit = true;
//
//const double convcrit = 1e-10;

//int main()
//{
//	// 555.50,867.97,1041.566
//	std::vector<double> vels = { 555.50, 867.97, 1041.566 };
//	std::vector<int> angles = { 30, 20, 10 };
//
//	for (auto vel : vels)
//	{
//		for (auto angle : angles)
//		{
//			std::clock_t c_start = std::clock();
//			executeSolver(vel, angle, eps, kappa, limit, convcrit);
//			std::clock_t c_end = std::clock();
//			double cpu_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
//			std::cout << "Run complete after " << cpu_time << "ms on CPU for v=" << std::to_string(vel) << " and alpha=" << std::to_string(angle) << "\n";
//		}
//	}
//
//	return 0;
//}




//int main()
//{
//	Solver s;
//	std::string meshpath = "../../../mesh/";
//	std::string respath = "../../../solution/";
//	std::string gridname = "sodX";
//	std::string infile = meshpath + gridname + ".grd";
//	std::string outfile = respath + gridname + ".res";
//
//	s.readGridGridPro(infile);
//	s.setSodXInitial();
//	int maxsteps = 1000;
//	for (int i = 1; i<=maxsteps;i++)
//	{
//		s.executeTimeStepGlobal();
//		if (i % 10 == 0)
//		{
//			outfile = respath + gridname + ".res";
//			outfile = respath + gridname + std::to_string(i) + ".res";
//			s.writeSolution(outfile);
//			std::cout << s.getResidual() << "\n";
//			//if (s.getResidual() < 1e-9)
//			//	break;
//		}
//	}
//	return 0;
//}
//
//int main()
//{
//	Solver s;
//	std::string meshpath = "../../../mesh/";
//	std::string respath = "../../../solution/";
//	std::string gridname = "Grid20deg";
//	std::string infile = meshpath + gridname + ".grd";
//	std::string outfile = respath + gridname + ".res";
//
//	s.readGridGridPro(infile);
//	// 555.50,867.97,1041.566
//	double vel = 555.50;
//	s.setConsInlet(1.01325e5, vel, 0, 300);
//	s.setConsInitial(1.01325e5, vel, 0, 300);
//	int maxsteps = 10000;
//	for (int i = 1; i <= maxsteps; i++)
//	{
//		s.executeTimeStepLocal();
//		if (i % (maxsteps / 1000) == 0)
//		{
//			outfile = respath + gridname + ".res";
//			outfile = respath + gridname + std::to_string(i) + ".res";
//			s.writeSolution(outfile);
//			std::cout << i << "\t" << s.getTime() << "\t" << s.getResiduals()[0] << "\n";
//			if (s.getResiduals()[0] < 1e-6)
//				break;
//		}
//	}
//	return 0;
//}