set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(boost_graph_ex tests.cpp exercices_toolbox.h)
    target_link_libraries(boost_graph_ex ${Boost_LIBRARIES} )
    ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
endif()
