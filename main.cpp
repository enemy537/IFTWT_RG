#include <iostream>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <exercices_toolbox.h>
#include <opencv2/opencv.hpp>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
<<<<<<< HEAD
#include <opencv2/imgproc/types_c.h>
=======
>>>>>>> a6f3865fa8202bf06b913435c1ce60d062cebfe9
#include "graph.hpp"

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointXYZI PointI;
typedef pcl::PointCloud<PointT> CloudT;
typedef pcl::PointCloud<PointI> CloudI;
<<<<<<< HEAD
=======
//typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS, pcl::PointXYZI, float> Graph;
>>>>>>> a6f3865fa8202bf06b913435c1ce60d062cebfe9

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
<<<<<<< HEAD
pcl::PointCloud<pcl::PointXYZI>::Ptr cloudGray(pcl::PointCloud<PointT>::Ptr cloud)
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
    return cloud_gray;
}

int main (int argc, char** argv) {
    CloudT::Ptr cloud_raw (new CloudT());
    if ( pcl::io::loadPCDFile <PointT> ("/home/pedro/Desktop/point_cloud/boamorte_final.pcd", *cloud_raw) == -1) {
=======
//float get_intensity_at(const pcl::PointCloud<PointI>::Ptr cloud, PointI& p){
//    pcl::KdTreeFLANN<PointI> tree;
//    tree.setInputCloud(cloud);
//    std::vector<int> nn_indices (1);
//    std::vector<float> nn_dists (1);
//
//    tree.nearestKSearch(p,1, nn_indices, nn_dists);
//
//    return cloud->points[nn_indices[0]].intensity;
//}
//auto static g_to_pc(Graph& g_)
//{
//    auto deep_copy = [](const PointI& p1){
//        PointI p (p1.intensity);
//        p.x = p1.x; p.y = p1.y; p.z = p1.z;
//        return p;
//    };
//
//    CloudI::Ptr out (new CloudI());
//    using vd = typename Graph::vertex_descriptor;
//    const auto vd_v = vertices(g_);
//    for (auto v = vd_v.first; v != vd_v.second; ++v) {
//
//        out->push_back(deep_copy(g_[*v]));
//    }
//    return out;
//};
//pcl::PointCloud<PointI>::Ptr erode_cloud(Graph& g, pcl::PointCloud<PointI>::Ptr c,
//                                         pcl::octree::OctreePointCloudAdjacency<PointI>& octree)
//{
//    using vd = typename Graph::vertex_descriptor;
//    pcl::PointCloud<PointI>::Ptr out (new pcl::PointCloud<PointI>);
//    Graph g_out;
//    octree.computeVoxelAdjacencyGraph(g_out);
//    using vd = typename Graph::vertex_descriptor;
//    const auto vd_g = vertices(g);
//    const auto vd_o = vertices(g_out);
//
//    int cnt =  0;
//    auto v = vd_g.first, v_o = vd_o.first;
//    for (; v != vd_g.second; ++v, ++v_o) {
//        std::cout << "Converting vertex " << ++cnt << std::endl;
//        const auto adj_v = boost::adjacent_vertices(*v, g);
//        bool erodible = std::find_if(adj_v.first, adj_v.second,
//                                     [&g](const vd &d) { return g[d].intensity == 0; }
//        ) != adj_v.second;
//        int intensity = erodible ? 0 : 255;
//        g_out[*v].intensity = intensity;
//    }
//
//    std::cout << "Nope " << std::endl;
//    return g_to_pc(g_out);
//}

int main (int argc, char** argv) {
    CloudT::Ptr cloud_raw (new CloudT());
    if ( pcl::io::loadPCDFile <PointT> ("/home/avell/Desktop/point_cloud/lucy.pcd", *cloud_raw) == -1) {
>>>>>>> a6f3865fa8202bf06b913435c1ce60d062cebfe9
        std::cout << "Cloud reading failed." << std::endl;
        return (-1);
    }

<<<<<<< HEAD
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_i = cloudGray(cloud_raw);

    Graph gg (cloud_i,10);

    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_e (new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_d (new pcl::PointCloud<pcl::PointXYZI>);
    gg.morph_erode(cloud_e);
    gg.morph_dilate(cloud_d);

    CloudI::Ptr cloud_g (new CloudI);
    pcl::copyPointCloud(*cloud_d,*cloud_g);

    for(int i = 0; i < cloud_g->size(); i++){
        float intensity = cloud_d->points[i].intensity
                          - cloud_e->points[i].intensity;
        cloud_g->points[i].intensity = intensity;
    }

//    pcl::io::savePCDFile("lucy_gradient.pcd",*cloud_g);
//    pcl::io::savePCDFile("lucy_erode.pcd",*cloud_e);
//    pcl::io::savePCDFile("lucy_dilate.pcd",*cloud_d);
=======
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_i = cloudRGB2GRAY(cloud_raw);

    pcl::VoxelGrid<PointI> grid;
    grid.setInputCloud(cloud_i);
    grid.setLeafSize(.006f,.006f,.006f);
    grid.filter(*cloud_i);

    Graph graph (cloud_i);

    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_d = graph.morph_bin_erode();
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_e = graph.morph_bin_dilate();

//    pcl::io::savePCDFile("lucy.pcd",*cloud_i);
//    pcl::io::savePCDFile("lucy_erode.pcd",*graph.morph_bin_erode());
//    pcl::io::savePCDFile("lucy_dilate.pcd",*graph.morph_bin_dilate());
>>>>>>> a6f3865fa8202bf06b913435c1ce60d062cebfe9

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

    // Normal visualization stuff
<<<<<<< HEAD
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3d cloud"));
    pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity (cloud_g,"intensity");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
    viewer->addPointCloud<pcl::PointXYZI>(cloud_g,intensity,"cloud");
=======
    /**
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3d cloud"));
    pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity (cloud_i,"intensity");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
    viewer->addPointCloud<pcl::PointXYZI>(cloud_i,intensity,"cloud");
    viewer->addCoordinateSystem (1.0);
>>>>>>> a6f3865fa8202bf06b913435c1ce60d062cebfe9
    viewer->setBackgroundColor (0, 0, 0);
    viewer->removeOrientationMarkerWidgetAxes();
    viewer->initCameraParameters ();

    while (!viewer->wasStopped ()){
        viewer->spinOnce(100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
<<<<<<< HEAD


=======
    **/
>>>>>>> a6f3865fa8202bf06b913435c1ce60d062cebfe9
    return (0);
}