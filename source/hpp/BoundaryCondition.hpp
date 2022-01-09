#pragma once
#include<Preprocessor.hpp>

enum BoundaryConditionType
{
	WALL = 1,
	INFLOW,
	OUTFLOW
};

class BoundaryCondition
{
public:
	BoundaryCondition();
	~BoundaryCondition();
	void SetBoundaryCondition(const int& type, const double& value = 0);
	int GetType();
	double GetValue();
private:
	int m_type = -1;
	double m_value = 0.0;
};
