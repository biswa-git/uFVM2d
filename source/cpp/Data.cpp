#pragma once
#include<Data.hpp>

EdgeData::EdgeData()
{
	m_edge_data.resize(EDGE_DATA_COUNT);
}

EdgeData::~EdgeData()
{
}

double& EdgeData::GetEdgeData(const int& data_id)
{
	return m_edge_data[data_id];
}

FaceData::FaceData()
{
	m_face_data.resize(FACE_DATA_COUNT);
}

FaceData::~FaceData()
{
}

double& FaceData::GetFaceData(const int& data_id)
{
	return m_face_data[data_id];
}





FaceGradData::FaceGradData()
{
	m_face_grad_data.resize(FACE_GRAD_DATA_COUNT);
	for (auto& face_grad_variable_data : m_face_grad_data)
	{
		face_grad_variable_data.resize(3);
	}
}

FaceGradData::~FaceGradData()
{
}

std::vector<double>& FaceGradData::GetFaceGradData(const int& data_id)
{
	return m_face_grad_data[data_id];
}