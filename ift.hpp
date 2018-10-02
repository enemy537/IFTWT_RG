#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
#include <queue>
#include <iostream>
#include <map>
#include <algorithm>

#define MAX_COST 500
#define MIN_COST 0


struct VertexBundle{
    unsigned int cost;
    PointI point;
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

typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS,
        VertexBundle,EdgeBundle> graph_v;
typedef boost::graph_traits<graph_v>::vertex_descriptor pcd_vx_descriptor;
typedef boost::graph_traits<graph_v>::vertex_iterator pcd_vx_iterator;


class IFT_PCD{
public:
    IFT_PCD(Graph* gg){
        g_t = gg->getAdjacencyList();
        copy_graph();
    }

    CloudI::Ptr gv_to_pc()
    {
        CloudI::Ptr out (new CloudI());
        const auto vd_v = vertices(g_v);
        for (auto v = vd_v.first; v != vd_v.second; ++v) {
            out->push_back(g_v[*v].point);
        }
        return out;
    };
private:

    void copy_graph(){
        std::map<graph_t::vertex_descriptor, pcd_vx_descriptor> map;
        BGL_FORALL_VERTICES(v,g_t,graph_t)
        {
            pcd_vx_descriptor pcd_vx = boost::add_vertex({MAX_COST,g_t[v]},g_v);
            map[v] = pcd_vx;
        }
        BGL_FORALL_EDGES(e,g_t,graph_t)
        {
            boost::add_edge(map[boost::source(e,g_t)],map[boost::target(e,g_t)]);
        }
    }

    graph_v g_v;
    graph_t g_t;
};