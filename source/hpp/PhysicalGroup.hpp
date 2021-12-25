#pragma once
#include<string>
#include<Face.hpp>

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
    void AddEdges(Edge*);
    void AddFaces(Face*);
private:
    PhysicalGroup(const int&, const std::string&);

    int m_type;
    std::string m_name;
    std::vector<Edge*> m_associated_edges;
    std::vector<Face*> m_associated_faces;
};
