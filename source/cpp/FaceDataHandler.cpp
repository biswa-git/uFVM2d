#include<FaceDataHandler.hpp>

FaceDataHandler::FaceDataHandler(Data& data):m_data(&data)
{
	m_face_list = data.GetFaceList();
}

FaceDataHandler::~FaceDataHandler()
{
}

void FaceDataHandler::UpdateGradient(const int& data_id)
{
	auto& edge_data = m_data->GetEdgeData(data_id);
	auto& face_grad_data = m_data->GetFaceGradData(data_id);

	for (auto& face : m_face_list)
	{
		auto face_id = face->GetId();
		auto& half_edges = face->GetHalfEdge();

		auto volume = face->GetArea();

		face_grad_data[face_id][0] = 0.0;
		face_grad_data[face_id][1] = 0.0;
		face_grad_data[face_id][2] = 0.0;

		for (auto half_edge : half_edges)
		{
			auto edge = half_edge->GetParentEdge();
			auto edge_id = edge->GetId();

			auto area_vector = half_edge->GetNormal()*half_edge->GetEdgeVector().Abs();

			face_grad_data[face_id][0] += edge_data[edge_id] * area_vector.GetDx();
			face_grad_data[face_id][1] += edge_data[edge_id] * area_vector.GetDy();
			face_grad_data[face_id][2] += edge_data[edge_id] * area_vector.GetDz();
		}
		face_grad_data[face_id][0] /= volume;
		face_grad_data[face_id][1] /= volume;
		face_grad_data[face_id][2] /= volume;
	}

}