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
#include <queue>
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

static bool p_equal (const PointI& pi, const Point& p){
    return pi.x == p.x && pi.y == p.y && pi.z == p.z;
};

class IFT_PCD{
public:
    IFT_PCD(graph_t& g_t, Cloud::Ptr seeds){
        this->g_t = g_t;
        this->seeds = seeds;
        copy_graph();
        add_seeds();
        compute_watershed();
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr getLabelCloud()
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr out (new pcl::PointCloud<pcl::PointXYZRGB>());

        std::map<pcd_vx_descriptor , std::vector<int>> colors;
        auto random = []{
            std::vector<int> colors;
            colors.push_back(rand()%255);
            colors.push_back(rand()%255);
            colors.push_back(rand()%255);
            return colors;
        };

        for(auto it : root_m){
            pcl::PointXYZRGB point;
            point.x = g_v[it.first].point.x;
            point.y = g_v[it.first].point.y;
            point.z = g_v[it.first].point.z;

            if(colors.find(it.second) == colors.end())
                colors[it.second] = random();

            std::vector<int> color = colors[it.second];
            point.r = color[0];
            point.g = color[1];
            point.b = color[2];

            out->points.emplace_back(point);
        }
        return out;
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
            boost::add_edge(map[boost::source(e,g_t)],map[boost::target(e,g_t)],g_v);
        }
    }
    void add_seeds()
    {
        graph_v g_in = g_v;
        BGL_FORALL_VERTICES(v,g_v,graph_v)
        {
            for(auto it = seeds->points.begin(); it != seeds->points.end(); it++)
            {
                bool pmin = p_equal(g_in[v].point,*it);
                if(pmin)
                {
                    g_in[v].cost = MIN_COST;
                    break;
                }
            }
        }
    }
    void compute_watershed()
    {
        /**
         * Trivial initialization
         */
        for(auto t : boost::make_iterator_range(boost::vertices(g_v)))
        {
            root_m[t] = t;
            if(g_v[t].cost != MAX_COST)
                Q.push(pair(t,g_v[t].cost));
        }

        /**
         * Propagation
         */
        while(!Q.empty())
        {
            pcd_vx_descriptor s = Q.top().first; Q.pop();
            const auto s_adj = boost::adjacent_vertices(s,g_v);
            for(auto t = s_adj.first; t != s_adj.second; ++t )
            {
                if(g_v[*t].cost > g_v[s].cost)
                {
                    float tmp = std::max(g_v[*t].point.intensity,(float)g_v[s].cost);
                    if(tmp < g_v[*t].cost)
                    {
                        if(g_v[*t].cost != MAX_COST)
                            Q.remove(pair(*t,g_v[*t].cost));
                        g_v[*t].cost = tmp;
                        root_m[*t] = root_m[s];
                        Q.push(pair(*t,g_v[*t].cost));
                    }
                }
            }
        }
    }

    graph_v g_v;
    graph_t g_t;
    Cloud::Ptr seeds;

    using pair = std::pair<pcd_vx_descriptor, unsigned int>;
    struct cmp {
        bool operator()(const pair &a, const pair &b) {
            return a.second > b.second;
        };
    };
    std::map<pcd_vx_descriptor , pcd_vx_descriptor> root_m;
    custom_priority_queue<pair, std::vector<pair>,cmp> Q;
};