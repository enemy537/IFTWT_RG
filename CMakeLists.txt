cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(IFTWT_RG)

set(CMAKE_CXX_STANDARD 11)


find_package(PCL 1.8 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


add_executable(IFTWT_RG main.cpp flat_point.hpp)
target_link_libraries (IFTWT_RG ${PCL_LIBRARIES})