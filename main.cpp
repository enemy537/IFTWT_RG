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
#include <opencv2/imgproc/types_c.h>
#include "graph.hpp"

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointXYZI PointI;
typedef pcl::PointCloud<PointT> CloudT;
typedef pcl::PointCloud<PointI> CloudI;

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
pcl::PointCloud<pcl::PointXYZI>::Ptr cloudGray(CloudT::Ptr cloud)
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
    if ( pcl::io::loadPCDFile <PointT> ("/home/pedro/Desktop/point_cloud/bmrt_small.pcd", *cloud_raw) == -1) {
        std::cout << "Cloud reading failed." << std::endl;
        return (-1);
    }

    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_i = cloudGray(cloud_raw);

    Graph gg (cloud_i,2);

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

    pcl::io::savePCDFile("bmrt_gradiente.pcd",*cloud_g);
//    pcl::io::savePCDFile("lucy_erode.pcd",*cloud_e);
//    pcl::io::savePCDFile("lucy_dilate.pcd",*cloud_d);

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
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3d cloud"));
    pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity (cloud_g,"intensity");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
    viewer->addPointCloud<pcl::PointXYZI>(cloud_g,intensity,"cloud");
    viewer->setBackgroundColor (0, 0, 0);
    viewer->removeOrientationMarkerWidgetAxes();
    viewer->initCameraParameters ();

    while (!viewer->wasStopped ()){
        viewer->spinOnce(100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
    return (0);
}