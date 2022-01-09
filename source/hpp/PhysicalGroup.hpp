#pragma once
#include<string>
#include<Face.hpp>
#include<BoundaryCondition.hpp>
enum PhysicalGroupType
{
    PG_BOUNDARY=1,
    PG_INTERIOR

};
class PhysicalGroup
{
public:
    ~PhysicalGroup();
    static PhysicalGroup* New(const int&, const std::string&);
    void SetBoundaryCondition(const BoundaryCondition& boundary_condition);
    BoundaryCondition PhysicalGroup::GetBoundaryCondition();
    std::string GetName();
    void AddEdges(Edge*);
    void AddFaces(Face*);
    std::vector<Edge*>& GetEdges();
private:
    PhysicalGroup(const int&, const std::string&);

    int m_type;
    BoundaryCondition m_boundary_condition;
    std::string m_name;
    std::vector<Edge*> m_associated_edges;
    std::vector<Face*> m_associated_faces;
};
