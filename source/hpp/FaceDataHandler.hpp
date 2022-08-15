#pragma once 
#include<Data.hpp>
#include<Face.hpp>

class FaceDataHandler
{
public:
	FaceDataHandler();
	~FaceDataHandler();
	void UpdateGradient(const int&);
private:
	std::vector<Face*> m_face_list;
};
