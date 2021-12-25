#include<Solution.hpp>

Solution::Solution()
{
	
}

Solution::~Solution()
{
}

GeometryResult Solution::ReadGeometry(const std::string& file_name)
{
	return geometry.Read(file_name);
}
