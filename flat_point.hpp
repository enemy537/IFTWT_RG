//w
// Created by avell on 16/03/18.
//
#include <pcl/segmentation/region_growing.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/extract_indices.h>

template <typename PointT, typename NormalT>
class FlatPoint : public pcl::RegionGrowing<PointT,NormalT>{
public:
    FlatPoint<PointT,NormalT>() :
        centroids_ (0),
        is_computed_centroid_info (false),
        centroids_cloud_ (new pcl::PointCloud<pcl::PointXYZ>),
        centroids_normals_ (new pcl::PointCloud<pcl::Normal>)
    {};

    std::string getClusterState(){
        return !this->clusters_.empty()?"Something":"Empty";
    }

    std::vector<Eigen::Vector4f> getCentroids() {
        if (this->centroids_.size() == 0){
            std::vector<pcl::PointIndices> clusters;

            if (this->clusters_.size() != 0)
                clusters = this->clusters_;
            else {
                //Compute the whole process
                this->extract(clusters);
            }

            for (const pcl::PointIndices &cluster : clusters) {
                Eigen::Vector4f centroid;
                pcl::compute3DCentroid(*this->input_, cluster, centroid);
                centroids_.push_back(centroid);
            }
        }
            // Vector of seeds
            return this->centroids_;
    }

    pcl::PointCloud<pcl::Normal>::Ptr getCentroidsNormals(){
        if(!is_computed_centroid_info)
            compute_centroid_info();
        return centroids_normals_;
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr getCentroidsCloud(){
        if(!is_computed_centroid_info)
            compute_centroid_info();
        return centroids_cloud_;
    }


protected:
    std::vector<Eigen::Vector4f> centroids_;
    boost::shared_ptr <pcl::PointCloud<pcl::PointXYZ> > centroids_cloud_;
    boost::shared_ptr <pcl::PointCloud<pcl::Normal> > centroids_normals_;
    bool is_computed_centroid_info;

    void compute_centroid_info(){
        getCentroids();

        std::vector<int> nn_indices (1);
        std::vector<float> nn_dists (1);

        for(const Eigen::Vector4f &centroid: centroids_ ){
            pcl::PointXYZ p;
            p.x = centroid[0]; p.y = centroid[1]; p.z = centroid[2];
            this->search_->nearestKSearch(p,1, nn_indices, nn_dists);

            centroids_cloud_->points.push_back(this->input_->points[nn_indices[0]]);

            nn_indices[0] = this->neighbour_number_;
            nn_dists[0]   = this->neighbour_number_;

            this->search_->nearestKSearch(p,this->neighbour_number_, nn_indices, nn_dists);

            for(const int &p_i: nn_indices){
                if(this->normals_->points[p_i].normal_x != 0) {
                    centroids_normals_->points.push_back(this->normals_->points[p_i]);
                    break;
                }
           }
        }
        is_computed_centroid_info = true;
    }
};

