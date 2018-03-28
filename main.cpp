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

int main (int argc, char** argv) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    if ( pcl::io::loadPCDFile <pcl::PointXYZ> ("/home/avell/Downloads/pcl-1.8.0/test/rops_cloud.pcd", *cloud) == -1) {
        std::cout << "Cloud reading failed." << std::endl;
        return (-1);
    }

    pcl::search::Search<pcl::PointXYZ>::Ptr tree = boost::shared_ptr<pcl::search::Search<pcl::PointXYZ> > (new pcl::search::KdTree<pcl::PointXYZ>);
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
    reg.setNumberOfNeighbours (10);
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

    return (0);
}