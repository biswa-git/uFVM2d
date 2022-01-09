#include<BoundaryCondition.hpp>

BoundaryCondition::BoundaryCondition()
{
}

BoundaryCondition::~BoundaryCondition()
{
}

void BoundaryCondition::SetBoundaryCondition(const int& type, const double& value)
{
	m_type = type;
	m_value = value;
}

int BoundaryCondition::GetType()
{
	return m_type;
}

double BoundaryCondition::GetValue()
{
	return m_value;
}