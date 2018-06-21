#include <iostream>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/extract_indices.h>
#include "flat_point.hpp"
#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <pcl/filters/voxel_grid.h>
#include <vtkRenderWindow.h>
#include <vtkCubeSource.h>
#include <vtkCleanPolyData.h>
#include <exercices_toolbox.h>
#include <opencv2/opencv.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <pcl/kdtree/flann.h>

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointCloud<PointT> CloudT;

pcl::PointCloud<pcl::PointXYZI>::Ptr cloudRGB2GRAY(pcl::PointCloud<PointT>::Ptr cloud)
{
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_gray (new pcl::PointCloud<pcl::PointXYZI>);
    cloud_gray->height = cloud->height;
    cloud_gray->width = cloud->width;

    for(pcl::PointCloud<PointT>::iterator it = cloud->begin(); it != cloud->end(); it++){
        // Color conversion
        cv::Mat pixel(1,1,CV_8UC3,cv::Scalar(it->r,it->g,it->b));
        cv::Mat temp;
        cv::cvtColor(pixel,temp,CV_RGB2GRAY);

        pcl::PointXYZI pointI;
        pointI.x = it->x;
        pointI.y = it->y;
        pointI.z = it->z;
        pointI.intensity = temp.at<uchar>(0,0);

        cloud_gray->push_back(pointI);
    }
    cv::Mat gray_values(1,cloud_gray->size(),CV_8U);
    cv::Mat temp;

    int counter = 0;
    for(pcl::PointCloud<pcl::PointXYZI>::iterator it = cloud_gray->begin(); it != cloud_gray->end(); it++){
        gray_values.at<uchar>(0,counter) = it->intensity;
        counter++;
    }
    double thres_v = cv::threshold(gray_values,temp,0,255,CV_THRESH_OTSU);
    std::cout << "Otsu threshold value = " << thres_v << std::endl;

    for(pcl::PointCloud<pcl::PointXYZI>::iterator it = cloud_gray->begin(); it != cloud_gray->end(); it++){
        float v = it->intensity;
        if(v < thres_v) {it->intensity = 0;  }
        else            {it->intensity = 255;}
    }

    return cloud_gray;
}
typedef pcl::PointXYZI Point;
float get_intensity_at(const pcl::PointCloud<Point>::Ptr cloud, Point& p){
    pcl::KdTreeFLANN<pcl::PointXYZI> tree;
    tree.setInputCloud(cloud);
    std::vector<int> nn_indices (1);
    std::vector<float> nn_dists (1);

    tree.nearestKSearch(p,1, nn_indices, nn_dists);

    return cloud->points[nn_indices[0]].intensity;;
}
typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS, pcl::PointXYZI, float> Graph;
typedef pcl::PointXYZI Point;
pcl::PointCloud<Point>::Ptr erode_cloud(Graph& g, pcl::PointCloud<Point>::Ptr c)
{
    using vd = typename Graph::vertex_descriptor;
    pcl::PointCloud<Point>::Ptr out (new pcl::PointCloud<Point>);
    auto deep_copy = [](const Point& p1){
        Point p (p1.intensity);
        p.x = p1.x; p.y = p1.y; p.z = p1.z;
        return p;
    };
    const auto vd_v = vertices(g);
    for(auto v = vd_v.first; v != vd_v.second; ++v)
    {
        const auto adj_v = boost::adjacent_vertices(*v,g);
        bool erodible = std::find_if(adj_v.first,adj_v.second,
                                 [&g,&c](const vd& d)
                                 {
                                 return get_intensity_at(c,g[d]) == 0;
                                 }
                                 )!= adj_v.second;
        int intensity = erodible ? 0:255;
        pcl::PointXYZI p = deep_copy(g[*v]);
        p.intensity = intensity;
        out->push_back(p);
    }
    return out;
}

int main (int argc, char** argv) {
    CloudT::Ptr cloud_raw (new CloudT());
    CloudT::Ptr cloud (new CloudT());
    if ( pcl::io::loadPCDFile <PointT> ("/home/avell/Desktop/point_cloud/lucy.pcd", *cloud_raw) == -1) {
        std::cout << "Cloud reading failed." << std::endl;
        return (-1);
    }

    // Cloud normalization

    pcl::VoxelGrid<PointT> downsampler;
    downsampler.setInputCloud(cloud_raw);
    downsampler.setLeafSize(.006f,.006f,.006f);
    downsampler.filter(*cloud);

    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_i = cloudRGB2GRAY(cloud);

    // Octree adjacency
    pcl::octree::OctreePointCloudAdjacency<pcl::PointXYZI> octree(.008f);
    octree.setInputCloud(cloud_i);
    octree.addPointsFromInputCloud();

    // Project IFTWT_RG
    /**pcl::search::Search<pcl::PointXYZ>::Ptr tree = boost::shared_ptr<pcl::search::Search<pcl::PointXYZ> > (new pcl::search::KdTree<pcl::PointXYZ>);
    pcl::PointCloud <pcl::Normal>::Ptr normals (new pcl::PointCloud <pcl::Normal>);
    pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> normal_estimator;
    normal_estimator.setSearchMethod (tree);
    normal_estimator.setInputCloud (cloud);
    normal_estimator.setKSearch (50);
    normal_estimator.compute (*normals);

    FlatPoint<pcl::PointXYZ, pcl::Normal> reg;
    reg.setMinClusterSize (50);
    reg.setMaxClusterSize (100000);
    reg.setSearchMethod (tree);
    reg.setNumberOfNeighbours (50);
    reg.setInputCloud (cloud);
    reg.setInputNormals (normals);
    reg.setSmoothnessThreshold (2.0 / 180.0 * M_PI);
    reg.setCurvatureThreshold (1.0);

    std::vector<Eigen::Vector4f> seeds = reg.getCentroids();

    std::vector <pcl::PointIndices> clusters;
    reg.extract (clusters);

    std::cout << "Number of points: " << cloud->points.size() << std::endl;
    std::cout << "Number of seeds: " << seeds.size() << std::endl;
    std::cout << "Number of clusters: " << clusters.size() << std::endl;

    for(const Eigen::Vector4f& s : seeds ){
        std::cout << s[0] <<" "<< s[1] <<" "<< s[2] << std::endl;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Visualization stuff
    pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud ();
    pcl::PointCloud <pcl::PointXYZ>::Ptr centroid_cloud = reg.getCentroidsCloud ();

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(colored_cloud);

    for(pcl::PointXYZ &point : centroid_cloud->points){
        viewer->addSphere(point,0.01,std::to_string(point.x));
    }
    viewer->addPointCloud<pcl::PointXYZRGB>(colored_cloud,rgb,"color_cloud");
    viewer->setBackgroundColor (255, 255, 255);
    viewer->setCameraPosition(-10,5,1,0,0,0,0);
    viewer->removeOrientationMarkerWidgetAxes();
    viewer->initCameraParameters ();


    while (!viewer->wasStopped ()){
        viewer->spinOnce(100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
    **/


    // Adjacency Graph
    typedef  boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS, pcl::PointXYZI, float> Graph;
    Graph adjacency_graph;
    octree.computeVoxelAdjacencyGraph(adjacency_graph);
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_2 = erode_cloud(adjacency_graph,cloud_i);

    /**
    Graph::vertex_iterator vertexIt, vertexEnd;
    Graph::adjacency_iterator neighbourIt, neighbourEnd;
    boost::tie(vertexIt, vertexEnd) = boost::vertices(adjacency_graph);
    int cnt = 0;
    for (; vertexIt != vertexEnd; ++vertexIt) {
        cout << *vertexIt << " is connected with ";
        boost::tie(neighbourIt, neighbourEnd) = boost::adjacent_vertices(*vertexIt, adjacency_graph);
        for (; neighbourIt != neighbourEnd; ++neighbourIt)
            cout << *neighbourIt << " ";
        cout << "\n"; cnt++;
    }
    cout << "Voxel: "<< cnt << endl;
    cout << "Points: "<< cloud->size() << endl;
    **/

    // Voxel visualization stuff
    /**
    typedef std::vector< PointT, Eigen::aligned_allocator<PointT> > AlignedPointTVector;
    AlignedPointTVector voxelCenters;

    voxelCenters.clear();

    octree.getOccupiedVoxelCenters(voxelCenters);
    double voxelSideLen = std::sqrt(octree.getVoxelSquaredSideLen());

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

    viewer->addPointCloud<PointT>(cloud, "color_cloud");
    viewer->setBackgroundColor (0, 0, 0);
    viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE,3,"color_cloud");
    viewer->setCameraPosition(0,0,0,1,0,1);

    while (!viewer->wasStopped ()){
        viewer->spinOnce();
    }
    **/

    // Normal visualization stuff
//    /**
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3d cloud"));
    pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity (cloud_2,"intensity");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
    viewer->addPointCloud<pcl::PointXYZI>(cloud_2,intensity,"cloud");
    viewer->addCoordinateSystem (1.0);
    viewer->setBackgroundColor (0, 0, 0);
    viewer->initCameraParameters ();

//    pcl::io::savePCDFile("lucy_bin.pcd",*cloud_i);
    pcl::io::savePCDFile("lucy_erode.pcd",*cloud_2);

    while (!viewer->wasStopped ()){
        viewer->spinOnce(100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
//    **/
    return (0);
}