#pragma once
#include<iostream>
#include<Solution.hpp>
#include<set>
#include<iostream>
#include<chrono>
#include<lis.h>
#include<PoissonSolver.hpp>

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
    /*
    solution.SetTimeStep(0.01);
    solution.SetDensity(1);
    solution.SetViscosity(0.01);
    solution.SetMaxTimeStep(10000);

    solution.SolveMomentum();

    solution.WriteSolution("test");
	*/

    //auto geometry = solution.GetGeometry();
	//solution.TestFeature();
	solution.Solve();
    /*PoissonSolver solver(geometry);
    auto& x = solver.GetResult();

    auto& faces = geometry.GetFaceList();
    double value;
    for (auto face : faces)
    {
        lis_vector_get_value(x, face->GetId(), &value);
        face->GetFaceData(F_TEMPERATURE) = value;
    }

	//------------------------------------------------------------------------
	auto vertex_list = geometry.GetVertexList();
	auto face_list = geometry.GetFaceList();
	std::ofstream myfile;
	auto file = "temperature1.dat";
	myfile.open(file);
	myfile << "TITLE = \"title\"\n";
	myfile << "VARIABLES = \"X\", \"Y\", \"Temperature\"\n";
	myfile << "ZONE N = " << vertex_list.size() << ", E = " << face_list.size() << ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE\n";
	myfile << "VARLOCATION = ([3] = CELLCENTERED)" << "\n";
	for (auto it : vertex_list)
	{
		auto coord = it->GetPositionVector();
		myfile << coord[0] << " ";
	}
	myfile << "\n";
	for (auto it : vertex_list)
	{
		auto coord = it->GetPositionVector();
		myfile << coord[1] << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		myfile << it->GetFaceData(F_TEMPERATURE) << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->GetHalfEdge()[0]->GetStart()->GetId() << " " << it->GetHalfEdge()[1]->GetStart()->GetId() << " " << it->GetHalfEdge()[2]->GetStart()->GetId() << "\n";
	}

	myfile.close();
	//------------------------------------------------------------------------
	*/
    return 0;

}