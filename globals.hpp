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

#include <ctime>
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

    template<class T,
            class Container = std::vector<T>,
            class Compare = std::less<typename Container::value_type>>
    class custom_priority_queue : public std::priority_queue<T,Container,Compare>
    {
    public:
        bool remove(const T& value)
        {
            auto it = std::find(this->c.begin(), this->c.end(), value);
            if(it != this->c.end())
            {
                this->c.erase(it);
                std::make_heap(this->c.begin(), this->c.end(), this->comp);
                return true;
            } else
                return false;
        }
        Container& getContainer() { return this->c; }
        void print()
        {
            for(auto it = this->c.begin(); it != this->c.end(); ++it)
                std::cout << *it << " ";
            std::cout << std::endl;
        }
    };

    static bool p_equal (const global::PointI p1, const global::PointI p2){
        return p1.x == p2.x && p1.y == p2.y && p1.z == p2.z;
    };
}

#endif //IFTWT_RG_GLOBALS_H
