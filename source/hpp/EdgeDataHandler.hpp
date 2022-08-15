#pragma once 
#include<Data.hpp>
#include<Edge.hpp>
class EdgeDataHandler
{
public:
	EdgeDataHandler();
	~EdgeDataHandler();
	void Update(const int&);
private:
	std::vector<Edge*> m_edge_list;
};
