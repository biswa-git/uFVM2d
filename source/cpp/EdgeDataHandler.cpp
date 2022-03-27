#include<EdgeDataHandler.hpp>

EdgeDataHandler::EdgeDataHandler(const Data& data): m_data(data)
{
	m_edge_list = data.GetEdgeList();
}

EdgeDataHandler::~EdgeDataHandler()
{
}

void EdgeDataHandler::Update(const int& data_id)
{
	auto& edge_data = m_data.GetEdgeData(data_id);
	auto& face_data = m_data.GetFaceData(data_id);

	std::vector<Face*> face(2, nullptr);

	int face_id[2];

	for (auto edge:m_edge_list)
	{
		int edge_id = edge->GetId();

		face[0] = edge->GetHalfEdge(0)->GetFace();
		face[1] = edge->GetHalfEdge(1)->GetFace();

		face_id[0] = face[0]->GetId();
		face_id[1] = face[1]->GetId();

		if (face[0] != nullptr && face[1] != nullptr)
		{
			edge_data[edge_id] = (face_data[face_id[0]] * face[1]->GetArea() + 
								  face_data[face_id[1]] * face[0]->GetArea()) /
								 (face[0]->GetArea() + face[1]->GetArea());
		}
	}
}
