#pragma once 
#include<Geometry.hpp>
#include<map>

enum SolutionData
{
	MASS_FLUX,
	MASS_FLUX_STAR,
	DATA_COUNT
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

private:
	std::vector<Edge*> m_edge_list;
	std::vector<Face*> m_face_list;
	std::vector<std::vector<double>> m_edge_data;
	std::vector<std::vector<double>> m_face_data;
};
