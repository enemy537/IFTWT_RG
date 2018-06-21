//
// Created by avell on 07/05/18.
//
#include <pcl/pcl_base.h>
#include <pcl/octree/octree.h>
#include <pcl/octree/octree_pointcloud_adjacency.h>
#include <pcl/segmentation/supervoxel_clustering.h>

template <typename PointT, typename NormalT>
class Graph : public pcl::Supervoxel<PointT>{
public:
    Graph<PointT>() :
        cloud_ (cloud)
    {};
protected:
    pcl::PointCloud<PointT>::Ptr cloud_;
private:

};