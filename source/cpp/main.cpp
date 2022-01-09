#include<iostream>
#include<Solution.hpp>
#include<set>
#include<iostream>
#include<chrono>
#include<lis.h>
int main()
{
    std::string name ="test";
    std::string file_location = "D:/Project/uFVM/input/";
    std::string file_name = name + ".msh";

    Solution solution;
    solution.ReadGeometry(file_location + file_name);

    solution.SetBoundaryCondition("inlet", INFLOW, 1.0);
	solution.SetBoundaryCondition("outlet", OUTFLOW);
	solution.SetBoundaryCondition("wall", WALL);
	solution.SetBoundaryCondition("cylinder", WALL);

    solution.SetTimeStep(0.01);
    solution.SetDensity(1);
    solution.SetViscosity(0.01);
    solution.SetMaxTimeStep(10000);

    solution.SolveMomentum();

    solution.WriteSolution("test");

    return 0;
}