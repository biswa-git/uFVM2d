cmake_minimum_required (VERSION 3.2)
project(solver LANGUAGES CXX)

find_package(OpenMP REQUIRED)

project (uFVM)
include_directories(${CMAKE_SOURCE_DIR}/source/hpp)
include_directories(${CMAKE_SOURCE_DIR}/lib/lis/include)

add_executable(FVM
source/cpp/main.cpp
source/cpp/Geometry.cpp
source/cpp/Vertex.cpp
source/cpp/Vector.cpp
source/cpp/Solution.cpp
source/cpp/Edge.cpp
source/cpp/Face.cpp
source/cpp/PhysicalGroup.cpp
source/cpp/BoundaryCondition.cpp
#source/cpp/EdgeDataHandler.cpp
#source/cpp/FaceDataHandler.cpp
source/cpp/Data.cpp
source/cpp/PoissonSolver.cpp
)

TARGET_LINK_LIBRARIES( FVM PUBLIC OpenMP::OpenMP_CXX ${CMAKE_SOURCE_DIR}/lib/lis/lis.lib)

install(TARGETS FVM DESTINATION ${CMAKE_SOURCE_DIR}/exe)
