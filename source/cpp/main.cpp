#pragma once
#include<iostream>
#include<Solution.hpp>
#include<set>
#include<iostream>
#include<chrono>
#include<lis.h>
#include<EdgeDataHandler.hpp>

int main()
{
    std::string name ="test";
    std::string file_location = "D:/Project/uFVM/input/";
    std::string file_name = name + ".msh";

    Solution solution;
    solution.ReadGeometry(file_location + file_name);
    /*
    solution.SetBoundaryCondition("inlet", INFLOW, 1.0);
	solution.SetBoundaryCondition("outlet", OUTFLOW);
	solution.SetBoundaryCondition("wall", WALL);
	solution.SetBoundaryCondition("cylinder", WALL);

    solution.SetTimeStep(0.01);
    solution.SetDensity(1);
    solution.SetViscosity(0.01);
    solution.SetMaxTimeStep(10000);

    solution.SolveMomentum();
    9681
    solution.WriteSolution("test");
    */

    auto geometry = solution.GetGeometry();
    Data data(geometry);
    auto& data1 = data.GetEdgeData(MASS_FLUX);
    data1[0] = 1000;
    EdgeDataHandler eh(data);
    eh.Update(MASS_FLUX);
    return 0;
}