//
// Created by pedro on 27/11/18.
//

struct mst_info{
    int region_1;
    int region_2;
    double weight;

    bool operator == (mst_info a){
        if(a.region_1 == region_1 && a.region_2 == region_2)
            return true;
        else
            return false;
    }
};
struct mst_cmp {
    bool operator()(const mst_info &a, const mst_info &b) {
        return a.weight > b.weight;
    };
};

class Cluster{
public:
    Cluster(std::map<std::pair<global::pcd_vx_descriptor,global::pcd_vx_descriptor>, double> mst, int num_regions){
        
    }

private:

    std::set<int> closest (){
        std::set<int> key;
        double value = INFINITY;
        for(auto& i : distance_map)
        {
            if(i.second < value)
            {
                key = i.first;
                value = i.second;
            }
        }
        return key;
    };
    bool key_compare (std::set<int> a, std::set<int> b){
        for(auto& a_ : a)
            for(auto& b_ : b)
                if(a_ == b_) return true;
        return false;
    };
    bool isKeyRoot (std::set<int> a, std::set<int> b, std::set<int> key ){
        std::set<int> root_1, root_2;
        std::set_difference(a.begin(),a.end(), key.begin(),key.end(),
                            std::inserter(root_1,root_1.begin()));
        std::set_difference(b.begin(),b.end(), key.begin(),key.end(),
                            std::inserter(root_2,root_2.begin()));

        return root_1 == root_2;
    }
    void update_key (std::set<int> a, std::set<int> b, std::set<int> key){
        double temp_weight = distance_map[a];
        std::set<int> temp_key (a);
        temp_key.insert(key.begin(),key.end());
        distance_map[temp_key] = temp_weight;
        distance_map.erase(b); distance_map.erase(a);
    }

    void print_key (std::set<int> key){
        for(auto& i : key)
            std::cout << i << " ";
        std::cout << "| " ; std::cout << distance_map[key] << std::endl;
    }


    std::map<std::set<int>, double> distance_map;
    std::vector<mst_info> TLC;
};

