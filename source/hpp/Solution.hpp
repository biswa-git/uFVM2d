#pragma once
#include<Geometry.hpp>
#include<Data.hpp>
#include <lis.h>
struct Parameter
{
	double dt = 0.0;
	double density = 1.0;
	double viscosity = 1.0/100.0;
	int max_timestep = 100;
	double gamma = 0.5;
};

class Solution
{
public:
	Solution();
	~Solution();
	GeometryResult ReadGeometry(const std::string& file_name);
	GeometryResult WriteSolution(const std::string& file_name);
	void SetBoundaryCondition(std::string name, const int& type, const double& value = 0);
	void SetTimeStep(const double& dt);
	void SetDensity(const double& density);
	void SetViscosity(const double& viscosity);
	void SetMaxTimeStep(const int& max_timestep);
	Geometry& GetGeometry();
	void SolveMomentum();
	void TestFeature(Data& data);
	void PoissonSolver();
private:
	Geometry geometry;
	Parameter parameter;
	int max_thread=8;

	LIS_MATRIX A_us, A_p, A_u;
	LIS_VECTOR b_us, b_vs, b_p, b_u, b_v, x;

	LIS_SOLVER solver_momentum, solver_pressure, solver_velocity;
};