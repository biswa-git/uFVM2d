#include<FaceDataHandler.hpp>

FaceDataHandler::FaceDataHandler()
{
}

FaceDataHandler::~FaceDataHandler()
{
}

void FaceDataHandler::UpdateGradient(const int& data_id)
{
	for (auto& face : m_face_list)
	{
		auto& half_edges = face->GetHalfEdge();
		auto volume = face->GetArea();

		face->GetFaceGradData(data_id)[0] = 0.0;
		face->GetFaceGradData(data_id)[1] = 0.0;
		face->GetFaceGradData(data_id)[2] = 0.0;

		for (auto half_edge : half_edges)
		{
			auto edge = half_edge->GetParentEdge();
			auto edge_id = edge->GetId();

			auto area_vector = half_edge->GetNormal()*half_edge->GetEdgeVector().Abs();

			face->GetFaceGradData(data_id)[0] += half_edge->GetParentEdge()->GetEdgeData(data_id) * area_vector.GetDx();
			face->GetFaceGradData(data_id)[1] += half_edge->GetParentEdge()->GetEdgeData(data_id) * area_vector.GetDy();
			face->GetFaceGradData(data_id)[2] += half_edge->GetParentEdge()->GetEdgeData(data_id) * area_vector.GetDz();
		}
		face->GetFaceGradData(data_id)[0] /= volume;
		face->GetFaceGradData(data_id)[1] /= volume;
		face->GetFaceGradData(data_id)[2] /= volume;
	}

}