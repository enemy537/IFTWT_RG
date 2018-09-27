#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
//#include "graph.hpp"
#include <queue>
#include <iostream>
#include <map>
#include <algorithm>

#define MAX_COST 500
#define MIN_COST 0


class VertexBundle{
public:
    VertexBundle( unsigned int cost = MAX_COST): cost(cost) { }

    friend bool operator > (const VertexBundle& v1, const VertexBundle& v2)
    {
        return v1.cost > v2.cost;
    }
    unsigned int cost;
    PointI* point;
};
class EdgeBundle {
public:
    EdgeBundle( unsigned int weight=0) : weight(weight){}

    unsigned int weight;
};

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
    void print()
    {
        for(auto it = this->c.begin(); it != this->c.end(); ++it)
            std::cout << *it << " ";
        std::cout << std::endl;
    }
};

struct vertex_interface{
    void operator()(PointI& v1, VertexBundle& v2){
        v2.point = &v1;
    }
};
struct edge_interface{
    template <typename Edge1, typename Edge2>
    void operator()(const Edge1& v1, const Edge2& v2){ }
};

class IFT_PCD{
    typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS,
            VertexBundle,EdgeBundle> graph_v;
    typedef boost::graph_traits<graph_v>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<graph_v>::vertex_iterator vertex_iterator;

public:
    IFT_PCD(Graph* gg){
        boost::copy_graph(gg,g,boost::vertex_copy(vertex_iterface())
        .edge_copy(edge_interface()));
    }

    CloudI::Ptr gv_to_pc()
    {
        auto deep_copy = [](const PointI& p1){
            PointI p (p1.intensity);
            p.x = p1.x; p.y = p1.y; p.z = p1.z;
            return p;
        };

        CloudI::Ptr out (new CloudI());
        const auto vd_v = vertices(g);
        for (auto v = vd_v.first; v != vd_v.second; ++v) {

            out->push_back(deep_copy(g[*v]));
        }
        return out;
    };
private:

    graph_v g;

};