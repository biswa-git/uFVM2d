#pragma once
#include <Vector.hpp>
#include <Vertex.hpp>
#include <Face.hpp>
#include <Data.hpp>

class HalfEdge;
class Face;
class Edge
{
public:
    static Edge* New(Vertex*, Vertex*);
    ~Edge();
    Edge(Edge const&) = delete;
    Edge& operator=(Edge const&) = delete;
    void Flip();
    size_t GetId();
    Vertex* GetStart();
    Vertex* GetEnd();
    HalfEdge* GetHalfEdge(const size_t&);
    Vector GetCenter();

    double& GetEdgeData(const int&);

    static size_t GetEdgeCount();
    static void Legalize(Edge*);
    static Vector DistanceVector(Vertex*, Vertex*);
    void SetAreaByDistance(const double&);
    double GetAreaByDistance();

private:
    Edge(Vertex*, Vertex*);

    static size_t m_count;
    size_t m_id;

    Vertex* m_start;
    Vertex* m_end;

    HalfEdge* m_half_edge[2];

    double m_area_by_distance;
    EdgeData m_edge_data;
};

class HalfEdge
{
public:
    static HalfEdge* New(Edge*, Vertex*);
    ~HalfEdge() = default;
    HalfEdge(HalfEdge const&) = delete;
    HalfEdge& operator=(HalfEdge const&) = delete;
    Edge* GetParentEdge();
    HalfEdge* GetNeighbourHalfEdge();
    Vertex* GetStart();
    Vertex* GetEnd();
    void SetNext(HalfEdge*);
    HalfEdge* GetNext();
    void SetFace(Face*);
    Face* GetFace();
    Vector GetEdgeVector();
    Vector GetNormal();
    Vector& GeAreaVector();
    void SetDirectionCoefficient(const int&);
    int GDirectionCoefficient();

private:
    HalfEdge(Edge*, Vertex*);
    Edge* m_parent;
    Vertex* m_start;
    Face* m_associated_Face;
    HalfEdge* m_next;
    Vector m_area_vector;
    int m_direction_coefficient;
};
