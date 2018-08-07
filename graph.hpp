//
// Created by avell on 07/05/18.
//
#include <pcl/pcl_base.h>
#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <pcl/search/kdtree.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <vtkRenderWindow.h>
#include <vtkCubeSource.h>
#include <vtkCleanPolyData.h>

typedef pcl::PointXYZI PointI;
typedef pcl::PointCloud<PointI> CloudI;
typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS, PointI, float> graph_t;

auto deep_copy = [](const PointI& p1){
    PointI p (p1.intensity);
    p.x = p1.x; p.y = p1.y; p.z = p1.z;
    return p;
};
auto static g_to_pc(graph_t& g_)
{
    CloudI::Ptr out (new CloudI());
    using vd = typename graph_t::vertex_descriptor;
    const auto vd_v = vertices(g_);
    for (auto v = vd_v.first; v != vd_v.second; ++v) {

        out->push_back(deep_copy(g_[*v]));
    }
    return out;
};
class Graph {
public:
    Graph(pcl::PointCloud<PointI>::Ptr cloud) : voxel_size(0),
                                                tree_(new pcl::search::KdTree<PointI>),
                                                octree_(nullptr)
    {
        cloud_ = cloud;
        tree_->setInputCloud(cloud_);
        initialize();
    };
    void setInputCloud(pcl::PointCloud<PointI >::Ptr cloud)
    {
        cloud_ = cloud;
        tree_->setInputCloud (cloud_);
    };
    graph_t getAdjacencyList()
    {
        return g_;
    };
    boost::shared_ptr<pcl::visualization::PCLVisualizer> voxelViewer()
    {
        //     Voxel visualization stuff
        typedef std::vector< PointI, Eigen::aligned_allocator<PointI> > AlignedPointTVector;
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

        viewer->addPointCloud<PointI>(cloud_, "color_cloud");
        viewer->setBackgroundColor (0, 0, 0);
        viewer->setPointCloudRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_POINT_SIZE,3,"color_cloud");
        return viewer;
    };
    CloudI::Ptr morph_bin_erode()
    {
        graph_t g_in = g_, g_out = copy_g();

        const auto vd_g = vertices(g_in);
        const auto vd_o = vertices(g_out);
#pragma omp parallel
        {
#pragma omp single
            {
                auto v = vd_g.first, v_o = vd_o.first;
                for (; v != vd_g.second; ++v, ++v_o) {
#pragma omp task
                    {
                        const auto adj_v = boost::adjacent_vertices(*v, g_in);
                        using vd = typename graph_t::vertex_descriptor;
                        bool erodible = std::find_if(adj_v.first, adj_v.second,
                                        [&g_in](const vd &d) { return g_in[d].intensity == 0; }
                        ) != adj_v.second;
                        int intensity = erodible ? 0 : 255;
                        g_out[*v].intensity = intensity;
                    }
                }
            }
        }

        return g_to_pc(g_out);
    }
    CloudI::Ptr morph_bin_dilate()
    {
        graph_t g_in = g_, g_out = copy_g();

        const auto vd_g = vertices(g_in);
        const auto vd_o = vertices(g_out);
#pragma omp parallel
        {
#pragma omp single
            {
                auto v = vd_g.first, v_o = vd_o.first;
                for (; v != vd_g.second; ++v, ++v_o) {
#pragma omp task
                    {
                        const auto adj_v = boost::adjacent_vertices(*v, g_in);
                        using vd = typename graph_t::vertex_descriptor;
                        bool erodible = std::find_if(adj_v.first, adj_v.second,
                                        [&g_in](const vd &d) { return g_in[d].intensity == 255; }
                        ) != adj_v.second;
                        int intensity = erodible ? 255 : 0;
                        g_out[*v].intensity = intensity;
                    }
                }
            }
        }

        return g_to_pc(g_out);
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
        }
        return res;
    }
    void optimizeGraph()
    {
        std::vector<int> nn_i (1);
        std::vector<float> nn_d (1);

        auto find_p = [tree = tree_,&nn_i,&nn_d](const PointI& p)
        {return tree->nearestKSearch(p,1, nn_i, nn_d);};

        const auto vd_v = vertices(g_);
        for(auto v = vd_v.first; v != vd_v.second; ++v)
        {
            const auto p = find_p(g_[*v]);
            g_[*v].intensity = cloud_->points[nn_i[0]].intensity;
        }
    }
private:
    pcl::PointCloud<PointI>::Ptr cloud_;
    pcl::search::KdTree<PointI>::Ptr tree_;
    double voxel_size;
    graph_t g_;
    pcl::octree::OctreePointCloudAdjacency<PointI>::Ptr octree_;

    void initialize()
    {
        voxel_size = computeCloudResolution();

        // Compute the octree adjacency tree
        octree_.reset(new pcl::octree::OctreePointCloudAdjacency<PointI>(voxel_size));
        octree_->setInputCloud(cloud_);
        octree_->addPointsFromInputCloud();

        std::cout << "octree_size " << octree_->size() << std::endl;
        std::cout << "cloud_size " << cloud_->size() << std::endl;

        // Adjacency Graph
        octree_->computeVoxelAdjacencyGraph(g_);
        optimizeGraph();
    };

    graph_t copy_g()
    {
        std::map<graph_t::vertex_descriptor, int> index;
        for (auto v : boost::make_iterator_range(boost::vertices(g_))) {
            index.insert(std::make_pair(v, index.size()));
        }
        graph_t g_out;
        boost::copy_graph(g_,g_out,
        boost::vertex_index_map(boost::make_assoc_property_map(index)));
        return g_out;
    }
};