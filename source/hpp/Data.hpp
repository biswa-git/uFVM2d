#pragma once 
#include<Geometry.hpp>
#include<map>

enum SolutionData
{
	MASS_FLUX,
	MASS_FLUX_STAR,
	U,
	U_S,
	V,
	V_S,
	TEMPERATURE,
	DATA_COUNT
};

enum FaceData
{
	U,
	U_S,
	V,
	V_S,
	TEMPERATURE,
	FACE_DATA_COUNT
};

enum EdgeData
{
	MASS_FLUX,
	MASS_FLUX_STAR,
	U,
	U_S,
	V,
	V_S,
	TEMPERATURE,
	EDGE_DATA_COUNT
};

enum FaceGradData
{
	PRESSURE,
	TEMPERATURE,
	GRAD_DATA_COUNT
};

class Data
{
public:
	Data(const Geometry&);
	~Data();
	const std::vector<Edge*>& GetEdgeList() const;
	const std::vector<Face*>& GetFaceList() const;
	std::vector<double>& GetEdgeData(const int&);
	std::vector<double>& GetFaceData(const int&);
	std::vector<std::vector<double>>& GetFaceGradData(const int&);

private:
	std::vector<Edge*> m_edge_list;
	std::vector<Face*> m_face_list;
	std::vector<std::vector<double>> m_edge_data;
	std::vector<std::vector<double>> m_face_data;
	std::vector<std::vector<std::vector<double>>> m_face_grad_data;
};
