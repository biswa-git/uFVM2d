#pragma once 
#include<Geometry.hpp>
#include<Data.hpp>
#include <lis.h>

class PoissonSolver
{
public:
	PoissonSolver();
	void SetSize(const LIS_INT& size);
	LIS_VECTOR& GetResult();
	void MatrixSetValue(const int& i, const int& j, const double& value);
	void BVectorSetValue(const int& i, const double& value);
	void XVectorSetValue(const int& i, const double& value);
	void Prepare();
	~PoissonSolver();

private:
	LIS_INT m_size;
	LIS_MATRIX m_A;
	LIS_VECTOR m_x;
	LIS_VECTOR m_b;
	LIS_SOLVER m_solver;
};
