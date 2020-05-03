//
// Created by pedro on 02/05/2020.
//

#ifndef IFTWT_RG_GRACO_HPP
#define IFTWT_RG_GRACO_HPP

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <pcl/search/kdtree.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <vtkRenderWindow.h>
#include <vtkCubeSource.h>
#include <vtkCleanPolyData.h>

enum class MORPH_COLOR : char {dilate = 3,erode = 4};

global::CloudT::Ptr gc_to_pc(global::graph_t& g_)
{
    auto deep_c_copy = [](const global::PointT& p1){
        global::PointT p;
        p.x = p1.x; p.y = p1.y; p.z = p1.z;
        p.r = p1.r; p.g = p1.g; p.b = p1.b;
        return p;
    };

    global::CloudT::Ptr out (new global::CloudT());
    using vd = typename global::graph_t::vertex_descriptor;
    const auto vd_v = vertices(g_);
    for (auto v = vd_v.first; v != vd_v.second; ++v) {

        out->push_back(deep_c_copy(g_[*v]));
    }
    return out;
};
class GraCo {
public:
    GraCo(global::CloudT::Ptr cloud) : voxel_size(0),overlap(1.01),
                                       tree_(new pcl::search::KdTree<global::PointT>),
                                       octree_(nullptr)
    {
        cloud_ = cloud;
        tree_->setInputCloud(cloud_);
        initialize();
    };
    GraCo(global::CloudT::Ptr cloud, double ovrlp) : voxel_size(0),
                                                     tree_(new pcl::search::KdTree<global::PointT>),
                                                     octree_(nullptr)
    {
        overlap = ovrlp;
        cloud_ = cloud;
        tree_->setInputCloud(cloud_);
        initialize();
    };
    ~GraCo(){ }
    void setInputCloud(global::CloudT::Ptr cloud)
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
        typedef std::vector< global::PointT, Eigen::aligned_allocator<global::PointT> > AlignedPointTVector;
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

        viewer->addPointCloud<global::PointT>(cloud_, "color_cloud");
        viewer->setBackgroundColor (0, 0, 0);
        viewer->setPointCloudRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_POINT_SIZE,3,"color_cloud");
        viewer->setCameraPosition(0,0,0,1,0,1);
        return viewer;
    };
    global::CloudT::Ptr morph_erode(global::CloudT::Ptr out)
    {
        return morph_c_base(MORPH_COLOR::erode, out);
    }
    global::CloudT::Ptr morph_dilate(global::CloudT::Ptr out)
    {
        return morph_c_base(MORPH_COLOR::dilate, out);
    }
//    global::CloudT::Ptr morph_gradient()
//    {
//        global::CloudT::Ptr cloud_e (new global::CloudT);
//        global::CloudT::Ptr cloud_d (new global::CloudT);
//        this->morph_erode(cloud_e);
//        this->morph_dilate(cloud_d);
//
//        global::CloudT::Ptr cloud_g (new global::CloudT);
//        pcl::copyPointCloud(*cloud_d,*cloud_g);
//
//        for(int i = 0; i < cloud_g->size(); i++){
//            float intensity = cloud_d->points[i].intensity
//                              - cloud_e->points[i].intensity;
//            cloud_g->points[i].intensity = intensity;
//        }
//        return cloud_g;
//    }
    global::graph_t pc_to_g(global::CloudT::Ptr cloud)
    {
        pcl::search::KdTree<global::PointT>::Ptr tree (new pcl::search::KdTree<global::PointT>);
        global::graph_t out = copy_g();
        tree->setInputCloud(cloud);

        BGL_FORALL_VERTICES(v,out,global::graph_t)
        {
            int id_out = find_p(tree,g_[v]);
            out[v].r = cloud->points[id_out].r;
            out[v].g = cloud->points[id_out].g;
            out[v].b = cloud->points[id_out].b;
        }
        return out;
    }
    global::CloudT::Ptr getCloud()
    {
        return gc_to_pc(g_);
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
                        int num_adj = 0, sum[3] = {0, 0, 0};
                        const auto adj_v = boost::adjacent_vertices(v, g_);
                        for (auto adj = adj_v.first; adj != adj_v.second; adj++) {
                            sum[0] = g_[*adj].r; sum[1] = g_[*adj].g; sum[2] = g_[*adj].b;
                            num_adj++;
                        }
                        g_[v].r = (int) sum[0] / num_adj;
                        g_[v].g = (int) sum[1] / num_adj;
                        g_[v].b = (int) sum[2] / num_adj;
                    }
                }
            }
        }
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
    int find_p(pcl::search::KdTree<global::PointT>::Ptr tree, const global::PointT& p)
    {
        std::vector<int> nn_i (1);
        std::vector<float> nn_d (1);
        tree->nearestKSearch(p,1, nn_i, nn_d);
        return nn_i[0];
    }
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
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
                        g_[*v].r = cloud_->points[nn_i].r;
                        g_[*v].g = cloud_->points[nn_i].g;
                        g_[*v].b = cloud_->points[nn_i].b;
                    }
                }
            }
        }
    }
#pragma clang diagnostic pop
private:
    pcl::PointCloud<global::PointT>::Ptr cloud_;
    pcl::search::KdTree<global::PointT>::Ptr tree_;
    double voxel_size, overlap;
    global::graph_t g_;
    pcl::octree::OctreePointCloudAdjacency<global::PointT>::Ptr octree_;

    void initialize()
    {
        voxel_size = computeCloudResolution();
        // Compute the octree adjacency tree
        octree_.reset(new pcl::octree::OctreePointCloudAdjacency<global::PointT>(voxel_size));
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
    global::CloudT::Ptr morph_c_base(MORPH_COLOR type, pcl::PointCloud<global::PointT>::Ptr &out) {
        global::graph_t g_in = g_;
        out->swap(*gc_to_pc(g_));
        pcl::search::KdTree<global::PointT>::Ptr tree (new pcl::search::KdTree<global::PointT>);
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
                        std::array<int, 3> morph_i;

                        if(type == MORPH_COLOR::erode){
                            morph_i = {255, 255, 255};
                            for(auto i = adj_v.first; i!=adj_v.second; i++){
                                int r = g_in[*i].r, g = g_in[*i].g, b = g_in[*i].b;
                                int temp_rgb = r+b+g, temp_i = morph_i[0]+morph_i[1]+morph_i[2];
                                if(temp_rgb < temp_i){
                                    morph_i[0] = r;
                                    morph_i[1] = g;
                                    morph_i[2] = b;
                                } ;
                            }
                        }else if(type == MORPH_COLOR::dilate){
                            morph_i = {0, 0, 0};
                            for(auto i = adj_v.first; i!=adj_v.second; i++){
                                int r = g_in[*i].r, g = g_in[*i].g, b = g_in[*i].b;
                                int temp_rgb = r+b+g, temp_i = morph_i[0]+morph_i[1]+morph_i[2];
                                if(temp_rgb > temp_i){
                                    morph_i[0] = r;
                                    morph_i[1] = g;
                                    morph_i[2] = b;
                                } ;
                            }
                        }

                        int id_out = find_p(tree,g_in[*v]);
                        out->points[id_out].r = morph_i[0];
                        out->points[id_out].g = morph_i[1];
                        out->points[id_out].b = morph_i[2];
                    }
                }
            }
        }
        return out;
    }
};

#pragma clang diagnostic pop

#endif //IFTWT_RG_GRACO_HPP
