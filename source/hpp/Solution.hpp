#pragma once
#include<Geometry.hpp>
#include<Data.hpp>
#include<PoissonSolver.hpp>
#include <lis.h>
#include<omp.h>
#include <unordered_map>

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
	void TestFeature();
	void UpdateEdgeData(const int& data_id);
	void UpdateFaceGradient(const int& data_id);
	void UpdateFlux(const int& vel_x_data_id, const int& vel_y_data_id, const int& flux_data_id);
	void CopyEdgeData(const int& source, const int& destination);
	void CopyFaceData(const int& source, const int& destination);
	void CopyFaceGradData(const int& source, const int& destination);
	void SetEdgeData(const int& data_id, const int& value);
	void SetFaceData(const int& data_id, const int& value);


	void Solve();
	void ConstructMomentumMatrix();
private:
	Geometry m_geometry;
	Parameter m_parameter;
	int m_max_thread = std::max(omp_get_max_threads()/2,1);

	std::unordered_map<int, int> m_edge_face_map;
	std::unordered_map<int, int> m_face_edge_map;
	std::unordered_map<int, int> m_face_grad_edge_map;

	PoissonSolver m_momentum;
	PoissonSolver m_pressure;
};