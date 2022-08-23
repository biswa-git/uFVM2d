#pragma once 
#include<map>
#include<vector>
enum FaceDataEnum
{
	F_VEL_X,
	F_VEL_Y,
	F_VEL_XS,
	F_VEL_YS,
	F_PRESSURE,
	F_PRESSURE_S,
	F_TEMPERATURE,
	F_CENTRAL_TERM,
	FACE_DATA_COUNT
};

enum EdgeDataEnum
{
	E_MASS_FLUX,
	E_MASS_FLUX_S,
	E_VEL_X,
	E_VEL_Y,
	E_VEL_XS,
	E_VEL_YS,
	E_TEMPERATURE,
	EDGE_DATA_COUNT
};

enum FaceGradDataEnum
{
	FG_TEMPERATURE,
	FACE_GRAD_DATA_COUNT
};

class EdgeData
{
public:
	EdgeData();
	~EdgeData();
	double& GetEdgeData(const int&);

private:
	std::vector<double> m_edge_data;
};

class FaceData
{
public:
	FaceData();
	~FaceData();
	double& GetFaceData(const int&);

private:
	std::vector<double> m_face_data;
};

class FaceGradData
{
public:
	FaceGradData();
	~FaceGradData();
	std::vector<double>& GetFaceGradData(const int&);

private:
	std::vector<std::vector<double>> m_face_grad_data;
};
