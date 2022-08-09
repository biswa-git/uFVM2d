#pragma once
#include<iostream>
#include<Solution.hpp>
#include<set>
#include<iostream>
#include<chrono>
#include<lis.h>
#include<EdgeDataHandler.hpp>
#include<PoissonSolver.hpp>

int main()
{
    std::string name ="test";
    std::string file_location = "D:/Project/uFVM/input/";
    std::string file_name = name + ".msh";

    Solution solution;
    solution.ReadGeometry(file_location + file_name);
    
    solution.SetBoundaryCondition("inlet", INFLOW, 250.0);
	solution.SetBoundaryCondition("outlet", INFLOW, 500.0);
	solution.SetBoundaryCondition("wall", INFLOW, 250.0);
	solution.SetBoundaryCondition("cylinder", INFLOW, 1000.0);
    /*
    solution.SetTimeStep(0.01);
    solution.SetDensity(1);
    solution.SetViscosity(0.01);
    solution.SetMaxTimeStep(10000);

    solution.SolveMomentum();

    solution.WriteSolution("test");
    */

    auto geometry = solution.GetGeometry();
    Data data(geometry);
    //auto& data1 = data.GetFaceData(MASS_FLUX);
    //EdgeDataHandler eh(data);
    //eh.Update(MASS_FLUX);
    //solution.TestFeature(data);

    PoissonSolver solver(geometry);
    auto& x = solver.GetResult();

    auto& faces = geometry.GetFaceList();
    auto face_temperature_data = data.GetFaceData(TEMPERATURE);
    double value;
    for (auto face : faces)
    {
        lis_vector_get_value(x, face->GetId(), &value);
        face_temperature_data[face->GetId()] = value;
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
		auto face_id = it->GetId();
		myfile << face_temperature_data[face_id] << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->GetHalfEdge()[0]->GetStart()->GetId() << " " << it->GetHalfEdge()[1]->GetStart()->GetId() << " " << it->GetHalfEdge()[2]->GetStart()->GetId() << "\n";
	}

	myfile.close();
	//------------------------------------------------------------------------

    return 0;

}