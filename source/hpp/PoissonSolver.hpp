#pragma once 
#include<Geometry.hpp>
#include<Data.hpp>
#include <lis.h>

class PoissonSolver
{
public:
	PoissonSolver(Geometry& geometry);
	LIS_VECTOR& GetResult();
	~PoissonSolver();

private:
	Geometry* m_geometry;
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR b;
	LIS_SOLVER solver;
};
