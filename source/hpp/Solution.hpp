#pragma once

#include<Geometry.hpp>

class Solution
{
public:
	Solution();
	~Solution();
	GeometryResult ReadGeometry(const std::string& file_name);

private:
	Geometry geometry;
};