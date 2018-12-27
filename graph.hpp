//
// Created by avell on 07/05/18.
//
#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <pcl/search/kdtree.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <vtkRenderWindow.h>
#include <vtkCubeSource.h>
#include <vtkCleanPolyData.h>

enum class MORPH : char {bin_erode = 1,bin_dilate = 2,dilate = 3,erode = 4};

global::CloudI::Ptr g_to_pc(global::graph_t& g_)
{
    auto deep_copy = [](const global::PointI& p1){
        global::PointI p (p1.intensity);
        p.x = p1.x; p.y = p1.y; p.z = p1.z;
        return p;
    };

    global::CloudI::Ptr out (new global::CloudI());
    using vd = typename global::graph_t::vertex_descriptor;
    const auto vd_v = vertices(g_);
    for (auto v = vd_v.first; v != vd_v.second; ++v) {

        out->push_back(deep_copy(g_[*v]));
    }
    return out;
};
class Graph {
public:
    Graph(global::CloudI::Ptr cloud) : voxel_size(0),overlap(1.01),
                                                tree_(new pcl::search::KdTree<global::PointI>),
                                                octree_(nullptr)
    {
        cloud_ = cloud;
        tree_->setInputCloud(cloud_);
        initialize();
    };
    Graph(global::CloudI::Ptr cloud, double ovrlp) : voxel_size(0),
                                                    tree_(new pcl::search::KdTree<global::PointI>),
                                                    octree_(nullptr)
    {
        overlap = ovrlp;
        cloud_ = cloud;
        tree_->setInputCloud(cloud_);
        initialize();
    };
    ~Graph(){ }
    void setInputCloud(global::CloudI::Ptr cloud)
    {
        cloud_ = cloud;
        tree_->setInputCloud (cloud_);
        initialize();
    };
    void setOverlap(double overlap)
    {
        this->overlap = overlap;
        initialize();
    }
    global::graph_t getAdjacencyList()
    {
        return g_;
    };
    boost::shared_ptr<pcl::visualization::PCLVisualizer> voxelViewer()
    {
        //     Voxel visualization stuff
        typedef std::vector< global::PointI, Eigen::aligned_allocator<global::PointI> > AlignedPointTVector;
        AlignedPointTVector voxelCenters;

        voxelCenters.clear();

        octree_->getOccupiedVoxelCenters(voxelCenters);
        double voxelSideLen = std::sqrt(octree_->getVoxelSquaredSideLen());

        boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));

        // Draw Voxels

        vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New ();

        // Create every cubes to be displayed

        for(int i = 0; i < voxelCenters.size (); i++) {
            auto s = static_cast<float> (voxelSideLen / 2.0);
            float x = voxelCenters[i].x;
            float y = voxelCenters[i].y;
            float z = voxelCenters[i].z;

            vtkSmartPointer<vtkCubeSource> wk_cubeSource = vtkSmartPointer<vtkCubeSource>::New();

            wk_cubeSource->SetBounds(x - s, x + s, y - s, y + s, z - s, z + s);
            wk_cubeSource->Update();

            appendFilter->AddInputData(wk_cubeSource->GetOutput());
        }
        // Remove any duplicate points
        vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New ();

        cleanFilter->SetInputConnection (appendFilter->GetOutputPort ());
        cleanFilter->Update ();

        //Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> multiMapper = vtkSmartPointer<vtkPolyDataMapper>::New ();

        multiMapper->SetInputConnection (cleanFilter->GetOutputPort ());

        vtkSmartPointer<vtkActor> multiActor = vtkSmartPointer<vtkActor>::New ();

        multiActor->SetMapper (multiMapper);

        multiActor->GetProperty ()->SetColor (1.0, 1.0, 1.0);
        multiActor->GetProperty ()->SetAmbient (1.0);
        multiActor->GetProperty ()->SetLineWidth (0.1);
        multiActor->GetProperty ()->EdgeVisibilityOn ();
        multiActor->GetProperty ()->SetOpacity (1.0);

        multiActor->GetProperty ()->SetRepresentationToWireframe ();

        viewer->getRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(multiActor);

        viewer->addPointCloud<global::PointI>(g_to_pc(g_), "color_cloud");
        viewer->setBackgroundColor (0, 0, 0);
        viewer->setPointCloudRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_POINT_SIZE,3,"color_cloud");
        viewer->setCameraPosition(0,0,0,1,0,1);
        return viewer;
    };
    global::CloudI::Ptr morph_bin_erode(global::CloudI::Ptr out)
    {
        return morph_base(MORPH::bin_erode,out);
    }
    global::CloudI::Ptr morph_bin_dilate(global::CloudI::Ptr out)
    {
        return morph_base(MORPH::bin_dilate,out);
    }
    global::CloudI::Ptr morph_erode(global::CloudI::Ptr out)
    {
        return morph_base(MORPH::erode,out);
    }
    global::CloudI::Ptr morph_dilate(global::CloudI::Ptr out)
    {
        return morph_base(MORPH::dilate,out);
    }
    global::CloudI::Ptr morph_gradient()
    {
        global::CloudI::Ptr cloud_e (new global::CloudI);
        global::CloudI::Ptr cloud_d (new global::CloudI);
        this->morph_erode(cloud_e);
        this->morph_dilate(cloud_d);

        global::CloudI::Ptr cloud_g (new global::CloudI);
        pcl::copyPointCloud(*cloud_d,*cloud_g);

        for(int i = 0; i < cloud_g->size(); i++){
            float intensity = cloud_d->points[i].intensity
                            - cloud_e->points[i].intensity;
            cloud_g->points[i].intensity = intensity;
        }
        return cloud_g;
    }
    global::graph_t pc_to_g(global::CloudI::Ptr cloud)
    {
        pcl::search::KdTree<global::PointI>::Ptr tree (new pcl::search::KdTree<global::PointI>);
        global::graph_t out = copy_g();
        tree->setInputCloud(cloud);

        BGL_FORALL_VERTICES(v,out,global::graph_t)
        {
            int id_out = find_p(tree,g_[v]);
            out[v].intensity = cloud->points[id_out].intensity;
        }
        return out;
    }
    global::CloudI::Ptr getCloud()
    {
        return g_to_pc(g_);
    }
    void smooth_by_mean()
    {
#pragma omp parallel
        {
#pragma omp single
            {
                BGL_FORALL_VERTICES(v, g_, global::graph_t)
                {
#pragma omp task
                    {
                        int num_adj = 0, sum = 0;
                        const auto adj_v = boost::adjacent_vertices(v, g_);
                        for (auto adj = adj_v.first; adj != adj_v.second; adj++) {
                            sum += g_[*adj].intensity;
                            num_adj++;
                        }
                        g_[v].intensity = (int) sum / num_adj;
                    }
                }
            }
        }
    }
    global::CloudI::Ptr morph_erode_gradient()
    {
        global::CloudI::Ptr out = morph_gradient();
        g_ = pc_to_g(out);
        return morph_erode(out);
    }
    global::Cloud::Ptr getRegmin()
    {
        global::Cloud::Ptr out(new global::Cloud);
#pragma omp parallel
        {
#pragma omp single
            {
                for (auto v : boost::make_iterator_range(boost::vertices(g_))) {
#pragma omp task shared(out)
                    {
                        float v_val = g_[v].intensity;
                        if (v_val == 0)
#pragma omp critical
                            out->push_back(global::pi2p(g_[v]));
                        else {
                            auto v_adj = boost::adjacent_vertices(v, g_);
                            std::vector<float> val_adj;
                            for (auto adj = v_adj.first; adj != v_adj.second; ++adj)
                                val_adj.emplace_back(g_[*adj].intensity);
                            if (v_val <= *std::min_element(val_adj.begin(), val_adj.end()))
#pragma omp critical
                                out->push_back(global::pi2p(g_[v]));
                        }
                    }
                }
            }
        }

        return out;
    }
protected:
    double
    computeCloudResolution ()
    {
        // Finding optimum voxel size
        double res = 0.0;
        int n_points = 0;
        int nres;
        std::vector<int> indices (2);
        std::vector<float> sqr_distances (2);

        for (size_t i = 0; i < cloud_->size (); ++i)
        {
            if (! pcl_isfinite ((*cloud_)[i].x))
            {
                continue;
            }
            //Considering the second neighbor since the first is the point itself.
            nres = tree_->nearestKSearch (i, 2, indices, sqr_distances);
            if (nres == 2)
            {
                res += sqrt (sqr_distances[1]);
                ++n_points;
            }
        }
        if (n_points != 0)
        {
            res /= n_points;
            // Overlap adjust
            res *= overlap;
        }
        return res;
    }
    int find_p(pcl::search::KdTree<global::PointI>::Ptr tree, const global::PointI& p)
    {
        std::vector<int> nn_i (1);
        std::vector<float> nn_d (1);
        tree->nearestKSearch(p,1, nn_i, nn_d);
        return nn_i[0];
    }
    void optimizeGraph() {
        const auto vd_v = vertices(g_);
#pragma omp parallel
        {
#pragma omp single
            {
                for (auto v = vd_v.first; v != vd_v.second; ++v) {
#pragma omp task
                    {
                        auto nn_i = find_p(tree_, g_[*v]);
                        g_[*v].intensity = cloud_->points[nn_i].intensity;
                    }
                }
            }
        }
    }
private:
    pcl::PointCloud<global::PointI>::Ptr cloud_;
    pcl::search::KdTree<global::PointI>::Ptr tree_;
    double voxel_size, overlap;
    global::graph_t g_;
    pcl::octree::OctreePointCloudAdjacency<global::PointI>::Ptr octree_;

    void initialize()
    {
        voxel_size = computeCloudResolution();
        // Compute the octree adjacency tree
        octree_.reset(new pcl::octree::OctreePointCloudAdjacency<global::PointI>(voxel_size));
        octree_->setInputCloud(cloud_);
        octree_->addPointsFromInputCloud();

        std::cout << "octree_size " << octree_->size() << std::endl;
        std::cout << "cloud_size " << cloud_->size() << std::endl;

        // Adjacency Graph
        octree_->computeVoxelAdjacencyGraph(g_);
        optimizeGraph();
        cloud_.reset();
    };
    global::graph_t copy_g()
    {
        std::map<global::graph_t::vertex_descriptor, int> index;
        for (auto v : boost::make_iterator_range(boost::vertices(g_))) {
            index.insert(std::make_pair(v, index.size()));
        }
        global::graph_t g_out;
        boost::copy_graph(g_,g_out,
                          boost::vertex_index_map(boost::make_assoc_property_map(index)));
        return g_out;
    }
    global::CloudI::Ptr morph_base(MORPH type, global::CloudI::Ptr& out){
        global::graph_t g_in = g_;
        out->swap(*g_to_pc(g_));
        pcl::search::KdTree<global::PointI>::Ptr tree (new pcl::search::KdTree<global::PointI>);
        tree->setInputCloud(out);
        const auto vd_g = vertices(g_in);
#pragma omp parallel
        {
#pragma omp single
            {
                for (auto v = vd_g.first; v != vd_g.second; ++v) {
#pragma omp task shared(out)
                    {
                        const auto adj_v = boost::adjacent_vertices(*v, g_in);
                        using vd = typename global::graph_t::vertex_descriptor;
                        float intensity = 0;

                        if(type == MORPH::bin_erode){
                            bool erodible = std::find_if(adj_v.first, adj_v.second,
                                                         [&g_in](const vd &d) { return g_in[d].intensity == 0; }
                            ) != adj_v.second;
                            intensity = erodible ? 0.f : 255.f;
                        }else if(type == MORPH::bin_dilate){
                            bool dilatable = std::find_if(adj_v.first, adj_v.second,
                                                          [&g_in](const vd &d) { return g_in[d].intensity == 255;}
                            ) != adj_v.second;
                            float intensity = dilatable ? 255.f:0.f;
                        }else if(type == MORPH::erode){
                            intensity = 255;
                            for(auto i = adj_v.first; i!=adj_v.second; i++){
                                float temp_i = g_in[*i].intensity;
                                if(temp_i < intensity) intensity = temp_i ;
                            }
                        }else if(type == MORPH::dilate){
                            intensity = 0;
                            for(auto i = adj_v.first; i!=adj_v.second; i++){
                                float temp_i = g_in[*i].intensity;
                                if(temp_i > intensity) intensity = temp_i ;
                            }
                        }

                        int id_out = find_p(tree,g_in[*v]);
                        out->points[id_out].intensity = intensity;
                    }
                }
            }
        }
        return out;
    }
};
