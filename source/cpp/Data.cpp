#pragma once
#include<Data.hpp>

Data::Data(const Geometry& geometry):m_edge_list(geometry.GetEdgeList()),m_face_list(geometry.GetFaceList())
{
	m_edge_data.resize(DATA_COUNT);
	m_face_data.resize(DATA_COUNT);
	m_face_grad_data.resize(DATA_COUNT);

	for (auto& edge_variable_data : m_edge_data)
	{
		edge_variable_data.resize(Edge::GetEdgeCount());
	}

	for (auto& face_variable_data : m_face_data)
	{
		face_variable_data.resize(m_face_list.size());
	}

	for (auto& face_grad_variable_data : m_face_grad_data)
	{
		face_grad_variable_data.resize(m_face_list.size());
		for (auto& data : face_grad_variable_data)
		{
			data.resize(3);
		}
	}
}

Data::~Data()
{
}

const std::vector<Edge*>& Data::GetEdgeList() const
{
	return m_edge_list;
}

const std::vector<Face*>& Data::GetFaceList() const
{
	return m_face_list;
}

std::vector<double>& Data::GetEdgeData(const int& data_id)
{
	return m_edge_data[data_id];
}

std::vector<double>& Data::GetFaceData(const int& data_id)
{
	return m_face_data[data_id];
}

std::vector<std::vector<double>>& Data::GetFaceGradData(const int& data_id)
{
	return m_face_grad_data[data_id];
}