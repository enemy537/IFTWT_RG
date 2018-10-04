#include <iostream>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <pcl/features/normal_3d_omp.h>
#include <exercices_toolbox.h>
#include <opencv2/opencv.hpp>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <opencv2/imgproc/types_c.h>
#include "graph.hpp"
#include "flat_point.hpp"
#include "ift.hpp"

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointCloud<PointT> CloudT;

CloudI::Ptr cloudRGB2GRAY(CloudT::Ptr cloud)
{
    CloudI::Ptr cloud_gray (new CloudI);
    cloud_gray->height = cloud->height;
    cloud_gray->width = cloud->width;

    for(CloudT::iterator it = cloud->begin(); it != cloud->end(); it++){
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
    for(CloudI::iterator it = cloud_gray->begin(); it != cloud_gray->end(); it++){
        gray_values.at<uchar>(0,counter) = it->intensity;
        counter++;
    }
    double thres_v = cv::threshold(gray_values,temp,0,255,CV_THRESH_OTSU);
    std::cout << "Otsu threshold value = " << thres_v << std::endl;

    for(CloudI::iterator it = cloud_gray->begin(); it != cloud_gray->end(); it++){
        float v = it->intensity;
        if(v < thres_v) {it->intensity = 0;  }
        else            {it->intensity = 255;}
    }

    return cloud_gray;
}
CloudI::Ptr cloudGray(CloudT::Ptr cloud)
{
    CloudI::Ptr cloud_gray (new CloudI);
    cloud_gray->height = cloud->height;
    cloud_gray->width = cloud->width;

    for(CloudT::iterator it = cloud->begin(); it != cloud->end(); it++){
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

    CloudI::Ptr cloud_i (new CloudI());
    if ( pcl::io::loadPLYFile <PointI> ("../cloud/igreja.ply", *cloud_i) == -1) {
        std::cout << "Cloud reading failed." << std::endl;
        return (-1);
    }

    Graph gg (cloud_i, 10);

    CloudI::Ptr cloud_g = gg.morph_gradient();
    graph_t g = gg.pc_to_g(cloud_g);

    pcl::search::Search<pcl::PointXYZI>::Ptr tree =
            boost::shared_ptr<pcl::search::Search<pcl::PointXYZI> > (new pcl::search::KdTree<pcl::PointXYZI>);
    pcl::PointCloud <pcl::Normal>::Ptr normals (new pcl::PointCloud <pcl::Normal>);
    pcl::NormalEstimationOMP<pcl::PointXYZI, pcl::Normal> normal_estimator;
    normal_estimator.setSearchMethod (tree);
    normal_estimator.setInputCloud (cloud_g);
    normal_estimator.setKSearch (50);
    normal_estimator.compute (*normals);

    FlatPoint<pcl::PointXYZI, pcl::Normal> reg;
    reg.setMinClusterSize (50);
    reg.setMaxClusterSize (100000);
    reg.setSearchMethod (tree);
    reg.setNumberOfNeighbours (50);
    reg.setInputCloud (cloud_g);
    reg.setInputNormals (normals);
    reg.setSmoothnessThreshold (5.0 / 180.0 * M_PI);
    reg.setCurvatureThreshold (1.0);

//    pcl::PointCloud <pcl::PointXYZ>::Ptr centroid_cloud = reg.getCentroidsCloud ();
//
//    IFT_PCD ift(g,centroid_cloud);
//
//    pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = ift.getLabelCloud();
//
//    pcl::io::savePLYFile("gradient.ply", *cloud_g);
//    pcl::io::savePLYFile("labels.ply", *colored_cloud);

// Visualization stuff

    std::vector <pcl::PointIndices> clusters;
    reg.extract (clusters);
    pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud ();

//    for(pcl::PointXYZ &point : centroid_cloud->points){
//        viewer->addSphere(point,.5,std::to_string(point.x));
//    }


    pcl::visualization::CloudViewer viewer ("Cluster viewer");
    viewer.showCloud(colored_cloud);
    while (!viewer.wasStopped ()) {}


    return (0);
}