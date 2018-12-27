#include "globals.hpp"
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <pcl/features/normal_3d_omp.h>
#include <opencv2/opencv.hpp>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <opencv2/imgproc/types_c.h>
#include "graph.hpp"
#include "flat_point.hpp"
#include "ift.hpp"

#include <cstdlib>
extern "C"{
#include "cluster.h" /* The C Clustering Library */
}

#include "hcluster.hpp"

global::CloudI::Ptr cloudGray(global::CloudT::Ptr cloud)
{
    global::CloudI::Ptr cloud_gray (new global::CloudI);
    cloud_gray->height = cloud->height;
    cloud_gray->width = cloud->width;

    for(global::CloudT::iterator it = cloud->begin(); it != cloud->end(); it++){
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

    global::CloudI::Ptr cloud_i (new global::CloudI());
    if ( pcl::io::loadPLYFile <global::PointI> ("../../clean/curia_sg.ply", *cloud_i) == -1) {
        std::cout << "Cloud reading failed main cloud." << std::endl;
        return (-1);
    }

    global::Cloud::Ptr cloud (new global::Cloud());
    if ( pcl::io::loadPLYFile <global::Point> ("../../rg/curia/seed_curia_4.ply", *cloud) == -1) {
        std::cout << "Cloud reading failed main cloud." << std::endl;
        return (-1);
    }

    Graph gg(cloud_i);

    IFT_PCD ift(gg,cloud);

    std::cout << "Now clustering..." << std::endl;

    Cluster c(ift.getMST(),ift.getRoots(),ift.getGraph(),14);

    global::CloudT::Ptr colored_cloud = c.getLabelCloud();

//    pcl::io::savePLYFile("curia_sg_min.ply", *colored_cloud);

    pcl::visualization::CloudViewer viewer("cloud");
    viewer.showCloud(colored_cloud);
    while (!viewer.wasStopped ()){}

    return (0);
}
