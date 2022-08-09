#pragma once 
#include<Data.hpp>

class FaceDataHandler
{
public:
	FaceDataHandler(Data&);
	~FaceDataHandler();
	void UpdateGradient(const int&);
private:
	std::vector<Face*> m_face_list;
	Data* m_data;
};
