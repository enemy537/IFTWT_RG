cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(IFTWT_RG)

set(CMAKE_CXX_STANDARD 11)


find_package(PCL 1.8 REQUIRED)
find_package(VTK REQUIRED)
find_package( OpenCV REQUIRED)

include(${VTK_USE_FILE})

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

include_directories(boost_helper)

add_executable(IFTWT_RG main.cpp flat_point.hpp graph.hpp)
target_link_libraries (IFTWT_RG ${PCL_LIBRARIES} ${VTK_LIBRARIES} ${OpenCV_LIBS})