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
void PhysicalGroup::AddEdgeFaceMap(const std::pair<Edge*, Face*>& edge_face_map)
{
    m_associated_edges_faces.emplace_back(edge_face_map);
}

void PhysicalGroup::AddFaces(Face* f)
{
    m_associated_faces.emplace_back(f);
}

std::vector<std::pair<Edge*, Face*>>& PhysicalGroup::GetEdgeFace()
{
    return m_associated_edges_faces;
}