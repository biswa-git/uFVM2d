#pragma once 
#include<Data.hpp>

class EdgeDataHandler
{
public:
	EdgeDataHandler(Data&);
	~EdgeDataHandler();
	void Update(const int&);
private:
	std::vector<Edge*> m_edge_list;
	Data* m_data;
};
