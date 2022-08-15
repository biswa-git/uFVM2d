#include<EdgeDataHandler.hpp>

EdgeDataHandler::EdgeDataHandler()
{
}

EdgeDataHandler::~EdgeDataHandler()
{
}

void EdgeDataHandler::Update(const int& data_id)
{
	std::vector<Face*> face(2, nullptr);

	for (auto edge:m_edge_list)
	{
		face[0] = edge->GetHalfEdge(0)->GetFace();
		face[1] = edge->GetHalfEdge(1)->GetFace();

		if (face[0] != nullptr && face[1] != nullptr)
		{
			edge->GetEdgeData(data_id) = (face[0]->GetFaceData(data_id) * face[1]->GetArea() +
								          face[1]->GetFaceData(data_id) * face[0]->GetArea()) /
								         (face[0]->GetArea() + face[1]->GetArea());
		}
	}
}
