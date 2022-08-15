#pragma once
#include<Edge.hpp>
#include <Data.hpp>
#include<vector>

class HalfEdge;

class Face
{
public:
    Face();
    virtual ~Face();
    Face(Face const&) = delete;
    Face& operator=(Face const&) = delete;
    virtual int GetId() = 0;
    virtual double GetArea() = 0;
    virtual Vector GetCentroid() = 0;
    virtual Vector GetNormalVector() = 0;
    virtual void SetOrphanedEdgeRemoveFlag(bool);
    virtual std::vector<HalfEdge*>& GetHalfEdge() = 0;
    virtual std::vector<Vector> GetVerticesVector() = 0;
    double& Face::GetFaceData(const int& data_id);
    std::vector<double>& Face::GetFaceGradData(const int& data_id);

protected:
    virtual void CalculateArea() = 0;
    static size_t m_count;
    size_t m_id;
    double m_area;
    Vector m_centroid;
    bool m_is_orphaned_edge_remove_flag;
    std::vector<HalfEdge*> m_half_edge;

    FaceData m_face_data;
    FaceGradData m_face_grad_data;
    

};

class TriFace :public Face
{
public:
    static TriFace* New(Vertex*, Vertex*, Vertex*, const Vector& normal, const size_t&);
    ~TriFace();
    virtual int GetId();
    virtual double GetArea();
    virtual Vector GetCentroid();
    virtual Vector GetNormalVector();
    virtual std::vector<HalfEdge*>& GetHalfEdge();
    virtual std::vector<Vector> GetVerticesVector();

protected:
    virtual void CalculateArea();
    virtual void CalculateCentroid();

private:
    TriFace(Vertex*, Vertex*, Vertex*, const Vector& normal,const size_t&);
};
