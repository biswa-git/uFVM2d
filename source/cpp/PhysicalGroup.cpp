#include<PhysicalGroup.hpp>


PhysicalGroup::PhysicalGroup(const int& type, const std::string& name) :m_type(type), m_name(name)
{};

PhysicalGroup::~PhysicalGroup()
{};

PhysicalGroup* PhysicalGroup::New(const int& type, const std::string& name)
{
    return new PhysicalGroup(type, name);
}

void PhysicalGroup::AddEdges(Edge* e)
{
    m_associated_edges.push_back(e);
}

void PhysicalGroup::AddFaces(Face* f)
{
    m_associated_faces.push_back(f);
}