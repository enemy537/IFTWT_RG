//
// Created by pedro on 27/11/18.
//

#ifndef IFTWT_RG_GLOBALS_H
#define IFTWT_RG_GLOBALS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
#include <pcl/pcl_base.h>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/pcd_io.h>

#include <math.h>
#include <set>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <queue>

namespace global{
    /** PCL types **/
    typedef pcl::PointXYZ Point;
    typedef pcl::PointCloud<Point> Cloud;
    typedef pcl::PointXYZRGB PointT;
    typedef pcl::PointCloud<PointT> CloudT;
    typedef pcl::PointXYZI PointI;
    typedef pcl::PointCloud<PointI> CloudI;
    typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS, PointI, float> graph_t;

    /** boost types **/
    struct VertexBundle{
        unsigned int cost;
        PointI point;
    };
    class EdgeBundle {
    public:
        EdgeBundle( unsigned int weight=0) : weight(weight){}

        unsigned int weight;
    };

    typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS,
    VertexBundle,EdgeBundle> graph_v;
    typedef boost::graph_traits<graph_v>::vertex_descriptor pcd_vx_descriptor;
    typedef boost::graph_traits<graph_v>::vertex_iterator pcd_vx_iterator;
}

#endif //IFTWT_RG_GLOBALS_H
