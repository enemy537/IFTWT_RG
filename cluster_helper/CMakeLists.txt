set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

add_library(
        cluster_helper
        cluster.h
        cluster.c
        data.h
        data.c
)
target_include_directories(cluster_helper
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(cluster_helper ${Boost_LIBRARIES} )
    ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
endif()
