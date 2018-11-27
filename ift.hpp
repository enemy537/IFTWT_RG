#define MAX_COST 500
#define MIN_COST 0

/** MST */
typedef struct {
    double height = 0.0;
    double area = 0.0;
    double volume = 0.0;
    bool status = false;

    void update(double seed_i){
        area += 1;
        if(seed_i > height)
            height = seed_i;
        volume += area*height;
    }

    void print(){
        std::cout << "height " << height << ", area " <<
                  area << ", volume " << volume << ", status " <<
                  status << std::endl;
    }
} r_info;

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

static bool p_equal (const global::PointI& pi, const global::Point& p){
    return pi.x == p.x && pi.y == p.y && pi.z == p.z;
};

class IFT_PCD{
public:
    IFT_PCD(global::graph_t& g_t, global::Cloud::Ptr seeds){
        this->g_t = g_t;
        this->seeds = seeds;
        copy_graph();
        add_seeds();
        compute_watershed();
    }

    global::CloudT::Ptr getLabelCloud()
    {
        global::CloudT::Ptr out (new global::CloudT());

        std::map<global::pcd_vx_descriptor , std::vector<int>> colors;
        auto random = []{
            std::vector<int> colors;
            colors.push_back(rand()%255);
            colors.push_back(rand()%255);
            colors.push_back(rand()%255);
            return colors;
        };

        for(auto it : root_m){
            global::PointT point;
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

    global::CloudI::Ptr gv_to_pc()
    {
        global::CloudI::Ptr out (new global::CloudI());
        const auto vd_v = vertices(g_v);
        for (auto v = vd_v.first; v != vd_v.second; ++v) {
            out->push_back(g_v[*v].point);
        }
        return out;
    };
private:

    void copy_graph(){
        std::map<global::graph_t::vertex_descriptor, global::pcd_vx_descriptor> map;
        BGL_FORALL_VERTICES(v,g_t,global::graph_t)
        {
            global::pcd_vx_descriptor pcd_vx = boost::add_vertex({MAX_COST,g_t[v]},g_v);
            map[v] = pcd_vx;
        }
        BGL_FORALL_EDGES(e,g_t,global::graph_t)
        {
            boost::add_edge(map[boost::source(e,g_t)],map[boost::target(e,g_t)],g_v);
        }
    }
    void add_seeds()
    {
        global::graph_v g_in = g_v;
        BGL_FORALL_VERTICES(v,g_v,global::graph_v)
        {
            for(auto it = seeds->points.begin(); it != seeds->points.end(); it++)
            {
                bool pmin = p_equal(g_in[v].point,*it);
                if(pmin)
                {
                    g_in[v].cost = MIN_COST;
                    /** MST */
                    double seed_i = g_in[v].point.intensity;
                    r_info_m[v] = r_info{seed_i,1};
                    break;
                }
            }
        }
    }

    void find_and_swap(global::pcd_vx_descriptor old_,  global::pcd_vx_descriptor new_)
    {
        for(auto it = root_mst.begin();it != root_mst.end();it++)
            if(it->second == old_)
                root_mst[it->first] = new_;
    }

    void compute_watershed()
    {
        /**
         * Trivial initialization
         */
        for(auto t : boost::make_iterator_range(boost::vertices(g_v)))
        {
            root_m[t] = t;
            root_mst[t] = t;
            if(g_v[t].cost != MAX_COST)
                Q.push(pair(t,g_v[t].cost));
        }

        /**
         * Propagation
         */
        while(!Q.empty())
        {
            global::pcd_vx_descriptor s = Q.top().first; Q.pop();
            /** MST */
            auto s_root = root_mst[s];
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

                        /** New member add - update region info **/
                        root_mst[*t] = root_mst[s];
                        r_info_m[s_root].update(g_v[*t].point.intensity);

                        const auto r_adj = boost::adjacent_vertices(*t,g_v);
                        for(auto r = r_adj.first; r != r_adj.second; ++r)
                        {
                            /** MST -----------BLOCK-START---------- **/
                            global::pcd_vx_descriptor r_root = root_mst[*r];
                            bool exists = false;
                            if(r_info_m.find(r_root) != r_info_m.end())
                                exists = true;

                            if(exists && r_root != s_root && !r_info_m[s_root].status)
                            {
                                r_info_m[s_root].status = true;
                                double s_volume = r_info_m[s_root].volume,
                                        r_volume = r_info_m[r_root].volume,
                                        MST_weight = 0;
                                if(s_volume < r_volume){
                                    MST_weight = s_volume;
                                    find_and_swap(s_root,r_root);

                                    r_info_m[r_root].area += r_info_m[s_root].area;
                                    r_info_m[r_root].volume += r_info_m[s_root].volume;
                                }else{
                                    MST_weight = r_volume;
                                    auto it = root_mst.find(r_root);
                                    find_and_swap(r_root,s_root);

                                    r_info_m[s_root].area += r_info_m[r_root].area;
                                    r_info_m[s_root].volume += r_info_m[r_root].volume;
                                }
                                mst[std::make_pair(s_root,r_root)] = MST_weight;
                            }
                            /** MST -----------BLOCK-END---------- **/
                        }
                    }
                }
            }
        }
    }

    global::graph_v g_v;
    global::graph_t g_t;
    global::Cloud::Ptr seeds;

    using pair = std::pair<global::pcd_vx_descriptor, unsigned int>;
    struct cmp {
        bool operator()(const pair &a, const pair &b) {
            return a.second > b.second;
        };
    };
    std::map<global::pcd_vx_descriptor, global::pcd_vx_descriptor> root_m;
    custom_priority_queue<pair, std::vector<pair>,cmp> Q;
    /** MST auxiliary structures **/
    std::map<std::pair<global::pcd_vx_descriptor,global::pcd_vx_descriptor>, double> mst;
    std::map<global::pcd_vx_descriptor, r_info> r_info_m;
    std::map<global::pcd_vx_descriptor, global::pcd_vx_descriptor> root_mst;

};