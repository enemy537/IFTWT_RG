//
// Created by avell on 18/06/18.
//

#ifndef BOOST_GRAPH_EX_EXERCICES_TOOLBOX_H
#define BOOST_GRAPH_EX_EXERCICES_TOOLBOX_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <fstream>

/** Trivial codes about vertices and edge manipulation
 *  without explicit values allocated
 * **/

boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::undirectedS>
create_empty_undirected_graph() noexcept
{
    return{};
};
boost::adjacency_list<>
create_empty_directed_graph() noexcept
{
    return{};
}
template <typename graph>
int get_n_vertices(const graph& g) noexcept
{
    const int n{
            static_cast<int>(boost::num_vertices(g))
    };
    assert(static_cast<unsigned long>(n) == boost::num_vertices(g));
    return n;
}
template <typename graph>
std::vector <
        typename boost::graph_traits<graph>::vertex_descriptor
            >
get_vertex_descriptors(const graph& g) noexcept
{
    using vd = typename graph::vertex_descriptor;

    std::vector<vd> vds(boost::num_vertices(g));
    const auto vis = vertices(g);
    std::copy(vis.first, vis.second,std::begin(vds));
    return vds;
}
template <typename graph>
std::pair<
        typename graph::vertex_iterator,
        typename graph::vertex_iterator
>
get_vertex_iterators(const graph& g) noexcept
{
    return vertices(g);
};
template <typename graph>
typename boost::graph_traits<graph>::vertex_descriptor
add_vertex(graph& g) noexcept
{
    static_assert(!std::is_const<graph>::value,
                  "graph cannot be const"
    );
    const auto vd = boost::add_vertex(g);
    return vd;
}
template <typename graph>
std::vector <
        typename boost::graph_traits<graph>::edge_descriptor
>
get_edge_descriptors(const graph& g) noexcept
{
    using ed = typename graph::edge_descriptor;

    std::vector<ed> v(boost::num_edges(g));
    const auto eip = edges(g);
    std::copy(eip.first, eip.second,std::begin(v));
    return v;
}
template <typename graph>
std::pair<
        typename graph::edge_iterator,
        typename graph::edge_iterator
>
get_edge_iterator(const graph& g) noexcept
{
    return edges(g);
};
template <typename graph>
int get_n_edges(const graph& g) noexcept
{
    const int n{
            static_cast<int>(boost::num_edges(g))
    };
    assert(static_cast<unsigned long>(n) == boost::num_edges(g));
    return n;
}
template <typename graph>
typename boost::graph_traits<graph>::edge_descriptor
add_edge(graph& g) noexcept
{
    static_assert(!std::is_const<graph>::value,
                  "graph cannot be const"
    );
    const auto vd_a = boost::add_vertex(g);
    const auto vd_b = boost::add_vertex(g);
    const auto aer = boost::add_edge(
            vd_a,vd_b,g
    );
    assert(aer.second);
    return aer.first;
}
boost::adjacency_list<>
create_markov_chain() noexcept
{
    auto g = create_empty_directed_graph();
    const auto vd_a = boost::add_vertex(g);
    const auto vd_b = boost::add_vertex(g);
    boost::add_edge(vd_a,vd_a,g);
    boost::add_edge(vd_b,vd_b,g);
    boost::add_edge(vd_a,vd_b,g);
    boost::add_edge(vd_b,vd_a,g);
    return g;
}
template <typename graph>
void save_graph_to_dot(graph& g, std::string& filename) noexcept
{
    std::ofstream f(filename);
    boost::write_graphviz(f,g);
}
boost::adjacency_list<boost::vecS, boost::vecS,
                      boost::undirectedS>
create_path_graph(const size_t n_vertices) noexcept
{
    auto g = create_empty_undirected_graph();
    if(n_vertices == 0) return g;
    auto vd_1 = boost::add_vertex(g);
    if(n_vertices == 1) return g;
    for(size_t i=1; i!=n_vertices; i++){
        auto vd_2 = boost::add_vertex(g);
        boost::add_edge(vd_1,vd_2,g);
        vd_1 = vd_2;
    }
    return g;
};
boost::adjacency_list<boost::vecS, boost::vecS,
                         boost::undirectedS>
create_k2_graph() noexcept
{
    auto g = create_empty_undirected_graph();
    const auto vd_a = boost::add_vertex(g);
    const auto vd_b = boost::add_vertex(g);
    boost::add_edge(vd_a,vd_b,g);
    return g;
};
boost::adjacency_list<boost::vecS, boost::vecS,
        boost::undirectedS>
create_peterson_graph() noexcept
{
    using vd = decltype(create_empty_undirected_graph())::vertex_descriptor;
    auto g = create_empty_undirected_graph();

    std::vector<vd> v;
    for (int i = 0; i!=5 ; ++i) {
        v.push_back(boost::add_vertex(g));
    }
    std::vector<vd> w;
    for (int i = 0; i!=5 ; ++i) {
        w.push_back(boost::add_vertex(g));
    }
    for (int i = 0; i!=5 ; ++i) {
        boost::add_edge(v[i],v[(1+i)%5],g);
    }
    for (int i = 0; i!=5 ; ++i){
        boost::add_edge(v[i],w[i],g);
    }
    for (int i = 0; i!=5 ; ++i){
        boost::add_edge(w[i],w[(i+2)%5],g);
    }
    return g;
};
template<typename graph>
bool has_edge_between_vertices(
        const typename boost::graph_traits<graph>::vertex_descriptor& vd_1,
        const typename boost::graph_traits<graph>::vertex_descriptor& vd_2,
        const graph& g
) noexcept
{
    return edge(vd_1,vd_2,g).second;
}
template<typename graph>
typename boost::graph_traits<graph>::edge_descriptor
get_edge_between_vertices(
        const typename boost::graph_traits<graph>::vertex_descriptor& vd_1,
        const typename boost::graph_traits<graph>::vertex_descriptor& vd_2,
        const graph& g
)
{
    const auto er = edge(vd_1,vd_2,g);
    if(!er.second)
    {
        std::stringstream msg;
        msg << __func__ << ": "
            << "no edge between these vertices";
        throw std::invalid_argument(msg.str());
    }
    return er.first;
}
template<typename graph, typename vertex_descriptor>
graph create_direct_neigbour_subgraph(
        const vertex_descriptor& vd,
        const graph& g
)
{
    graph h;
    std::map<vertex_descriptor,vertex_descriptor> m;
    {
        const auto vd_h = boost::add_vertex(h);
        m.insert(std::make_pair(vd,vd_h));
    }
    // Copy vertices
    {
        const auto vdsi = boost::adjacent_vertices(vd,g);
        for(auto i = vdsi.first; i != vdsi.second; i++)
        {
            if(m.find(*i) == m.end())
            {
                const auto vd_h = boost::add_vertex(h);
                m.insert(std::make_pair(*i,vd_h));
            }
        }
    }
    // Copy edges
    {
        const auto eip = edges(g);
        const auto j = eip.second;
        for (auto i=eip.first; i!=j; ++i)
        {
            const auto vd_from = source(*i,g);
            const auto vd_to = target(*i,g);
            if(m.find(vd_from) == std::end(m)) continue;
            if(m.find(vd_to) == std::end(m)) continue;
            boost::add_edge(m[vd_from],m[vd_to],h);
        }
    }
    return h;
};

/** Working with bundled values!
 *  This concept is applied to a sample class
 *  but can be extented to any other
 * **/
struct my_bundled_vertex{
    explicit my_bundled_vertex(
        const std::string& name = "",
        const std::string& description = "",
        const double x = 0.0,
        const double y = 0.0
    ) noexcept ;
    std::string m_name;
    std::string m_description;
    double m_x;
    double m_y;
};
std::ostream& operator << (std::ostream& os, const my_bundled_vertex& e) noexcept;
bool operator == (const my_bundled_vertex& lhs, const my_bundled_vertex& rhs) noexcept;
bool operator != (const my_bundled_vertex& lhs, const my_bundled_vertex& rhs) noexcept;

boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::directedS,
    my_bundled_vertex
>
create_empty_directed_bundled_vertices_graph() noexcept
{
    return {};
};
boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::undirectedS,
        my_bundled_vertex
>
create_empty_undirected_bundled_vertices_graph() noexcept
{
    return {};
};
template <typename graph, typename bundled_vertex>
typename boost::graph_traits<graph>::vertex_descriptor
add_bundled_vertex(const bundled_vertex& v, graph& g) noexcept
{
    static_assert(!std::is_const<graph>::value,
                  "graph cannot be const"
    );
    return boost::add_vertex(v,g);
}
template <typename graph>
std::vector<pcl::PointXYZI> get_my_bundled_vertexes(
        const graph& g) noexcept
{
    using vd = typename graph::vertex_descriptor;
    std::vector<pcl::PointXYZI> v(boost::num_vertices(g));
    const auto vip = vertices(g);
    std::transform(vip.first,vip.second,std::begin(v),
                   [&g](const vd& d){ return g[d]; }
    );
    return v;
};
#endif //BOOST_GRAPH_EX_EXERCICES_TOOLBOX_H
