cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(IFTWT_RG)

set(CMAKE_CXX_STANDARD 14)

find_package(PCL 1.8 REQUIRED)
find_package(VTK REQUIRED)
find_package( OpenCV REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
if (APPLE)
    FIND_PATH( GLEW_INCLUDE_DIR glew.h
            /System/Library/Frameworks/GLEW.framework/Versions/A/Headers
            ${OPENGL_LIBRARY_DIR}
            )
    SET(GLEW_GLEW_LIBRARY "-lGLEW" CACHE STRING "GLEW library for OSX")
    SET(GLEW_cocoa_LIBRARY "-framework Cocoa" CACHE STRING "Cocoa framework for OSX")
endif (APPLE)

include(${VTK_USE_FILE})

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

include_directories(boost_helper)

add_executable(IFTWT_RG main.cpp flat_point.hpp graph.hpp ift.hpp)
target_link_libraries (IFTWT_RG ${PCL_LIBRARIES} ${VTK_LIBRARIES} ${OpenCV_LIBS})