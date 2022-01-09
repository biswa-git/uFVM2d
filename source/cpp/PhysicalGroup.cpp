#include<PhysicalGroup.hpp>


PhysicalGroup::PhysicalGroup(const int& type, const std::string& name) :m_type(type), m_name(name)
{};

PhysicalGroup::~PhysicalGroup()
{};

PhysicalGroup* PhysicalGroup::New(const int& type, const std::string& name)
{
    return new PhysicalGroup(type, name);
}

void PhysicalGroup::SetBoundaryCondition(const BoundaryCondition& boundary_condition)
{
    m_boundary_condition = boundary_condition;
}

BoundaryCondition PhysicalGroup::GetBoundaryCondition()
{
    return m_boundary_condition;
}

std::string PhysicalGroup::GetName()
{
    return m_name;
}
void PhysicalGroup::AddEdges(Edge* e)
{
    m_associated_edges.push_back(e);
}

void PhysicalGroup::AddFaces(Face* f)
{
    m_associated_faces.push_back(f);
}

std::vector<Edge*>& PhysicalGroup::GetEdges()
{
    return m_associated_edges;
}