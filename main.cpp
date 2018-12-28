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
void
setBackground (pcl::visualization::PCLVisualizer& viewer)
{
    viewer.setBackgroundColor (.0, .0, .0);
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

    Graph gg(cloud_i,6);

    IFT_PCD ift(gg,cloud);

    std::cout << "Now clustering..." << std::endl;

    Cluster c(ift.getMST(),ift.getRoots(),ift.getGraph(),14);

    global::CloudT::Ptr colored_cloud_tree = c.getLabelCloud();
    global::CloudT::Ptr colored_cloud_hash = c.getCloudFromHash(14);

//    pcl::io::savePLYFile("curia_sg_min.ply", *colored_cloud);

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->initCameraParameters ();

    int v1(0);
    viewer->createViewPort (0.0, 0.0, 0.5, 1.0, v1);
    viewer->setBackgroundColor (0, 0, 0, v1);
    viewer->addText ("tree", 10, 10, "v1 text", v1);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb (colored_cloud_tree);
    viewer->addPointCloud<pcl::PointXYZRGB> (colored_cloud_tree, rgb, "sample cloud1", v1);

    int v2(0);
    viewer->createViewPort (0.5, 0.0, 1.0, 1.0, v2);
    viewer->setBackgroundColor (0.3, 0.3, 0.3, v2);
    viewer->addText ("hash", 10, 10, "v2 text", v2);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb2 (colored_cloud_hash);
    viewer->addPointCloud<pcl::PointXYZRGB> (colored_cloud_hash, rgb2, "sample cloud2", v2);

    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "sample cloud1");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "sample cloud2");
    viewer->addCoordinateSystem (1.0);

    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }

    return (0);
}
