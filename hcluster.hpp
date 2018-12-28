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
    Cluster(std::map<std::pair<global::pcd_vx_descriptor,global::pcd_vx_descriptor>, double> mst,
            std::map<global::pcd_vx_descriptor, global::pcd_vx_descriptor> root_m,
            global::graph_v g_v, int num_regions){

        this->g_v = g_v;
        this->root_m = root_m;
        create_TLC(mst);
        std::cout << "TLC = " << TLC.size() << std::endl;
        clusterize(num_regions);
    }

    global::CloudT::Ptr getLabelCloud()
    {
        global::CloudT::Ptr out (new global::CloudT());

        std::cout << "regions: " << regions.size() << std::endl;
        std::cout << "root_m: " << root_m.size() << std::endl;

        std::map<std::set<int>, std::vector<int>> colors;
        auto random = [] {
            std::vector<int> colors;
            colors.emplace_back(rand() % 255);
            colors.emplace_back(rand() % 255);
            colors.emplace_back(rand() % 255);
            return colors;
        };

        for(auto& it : regions) colors[it] = random();
        for(auto it : root_m){
            global::PointT point;
            point.x = g_v[it.first].point.x;
            point.y = g_v[it.first].point.y;
            point.z = g_v[it.first].point.z;
            std::set<int> key;
            int tlc_idx;

            int r = root_translator[it.second];
            for(int i = 0; i < TLC.size(); i++){
                if(r == TLC[i].region_1 || r == TLC[i].region_2){
                    tlc_idx = i; break;
                }
            }

            for(auto& c : regions){
                if(c.find(tlc_idx)!=c.end()){
                    key = c;
                    break;
                }
            }
            std::vector<int> color = colors[key];
            point.r = color[0]; point.g = color[1];point.b = color[2];

            out->points.emplace_back(point);
        }
        return out;
    }
    global::CloudT::Ptr getCloudFromHash(int num_regions){
        compute_distance_map();
        regions.clear();
        single_linking(num_regions);
        return getLabelCloud();
    }
private:
    void create_TLC(std::map<std::pair<global::pcd_vx_descriptor,global::pcd_vx_descriptor>, double> mst){

        int counter = 0;
        for (const auto& item: root_m) {
            if (root_translator.find(item.second) == root_translator.end())
                root_translator[item.second] = counter++;
        }
        global::custom_priority_queue<mst_info, std::vector<mst_info>,mst_cmp> mst_translator;

        for (const auto& item: mst)
            mst_translator.push(mst_info{root_translator[item.first.first],
                                         root_translator[item.first.second],
                                         item.second});
        while(!mst_translator.empty())
        {
            mst_info temp = mst_translator.top(); mst_translator.pop();
            TLC.push_back(temp);
            std::vector<mst_info> container = mst_translator.getContainer();
            for(auto it = container.begin(); it != container.end(); it++)
            {
                int r1 = it->region_1, r2 = it->region_2;
                if(r1 == temp.region_1){it->region_1 = counter;}
                else if(r1 == temp.region_2){it->region_2 = counter;}
                if(r2 == temp.region_1){it->region_1 = counter;}
                else if(r2 == temp.region_2){it->region_2 = counter;}
            }
            counter++;
        }
    }

    void clusterize(int num_regions){
        int nrows = TLC.size(), ncols = 1;
        double** data = (double**) malloc(nrows*sizeof(double*));
        int** mask = (int**) malloc(nrows*sizeof(int*));
        double* weight = (double*) malloc(ncols*sizeof(double));

        for (int i = 0; i < ncols; i++) weight[i] = 1.0;
        for (int i = 0; i < nrows; i++){
            data[i] = (double*) malloc(sizeof(double));
            mask[i] = (int*) malloc(sizeof(int));
            data[i][0] = TLC[i].weight;
            mask[i][0] = 1;
        };

        double** matrix = (double**) malloc(nrows*sizeof(double*));
        if (matrix) {
            matrix[0] = NULL;
            for (int i = 1; i < nrows; i++)
            { matrix[i] = (double*) malloc(i*sizeof(double));
                if (matrix[i]==NULL) /* Not enough memory available */
                { while (--i >= 0) free(matrix[i]);
                    free(matrix);
                    matrix = NULL;
                    break;
                }
            }
        }
        if (!matrix) {
            free(weight);
            printf("Insufficient memory to store the distance matrix\n");
        }

        distancematrix(nrows, ncols, data, mask, weight, 'e', 0, matrix);

        Node* tree = treecluster(nrows, ncols, 0, 0, 0, 0, 'e', 's', matrix);

        int* clusterid = (int*) malloc(nrows*sizeof(int));
        cuttree(nrows, tree, num_regions, clusterid);

        for(int i=0; i<num_regions; i++)
            regions.emplace_back(std::set<int>());
        for(int i=0; i<nrows; i++)
            regions[clusterid[i]].insert(i);
    }
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

    void compute_distance_map(){
        for (int i = 0; i < TLC.size(); i++) {
            for (int j = 0; j < TLC.size(); j++) {
                if (i != j && i < j)
                    distance_map[{i, j}] = std::sqrt(std::pow(TLC[i].weight - TLC[j].weight, 2));
            }
        }
    }

    void single_linking(int num_regions){
        for(int i = 0; i < (TLC.size() - num_regions); i++)
        {
            std::set<int> key = closest();
            std::vector<std::set<int>> keys;
            for(auto& d : distance_map)
                if(key_compare(key,d.first) && d.first != key){
                    std::set<int> deep_copy(d.first);
                    keys.emplace_back(deep_copy);
                }

            auto it_1 = keys.begin();
            auto it_2 = keys.begin();

            while(it_1 != keys.end()){
                bool deleted = false;
                while(it_2 != keys.end()){
                    if(isKeyRoot(*it_1,*it_2,key) && *it_1 != *it_2)
                    {
                        if(distance_map[*it_1] < distance_map[*it_2])
                            update_key(*it_1,*it_2,key);
                        else
                            update_key(*it_2,*it_1,key);

                        keys.erase(it_2); keys.erase(it_1);
                        it_1 = keys.begin(); it_2 = keys.begin();
                        deleted = true;
                        break;
                    }else it_2++;
                }
                if(!deleted)
                    it_1++;
            }
            distance_map.erase(key);
        }

        regions.emplace_back(distance_map.begin()->first);
        for(auto& i : distance_map){
            std::set<int> r1, r2, r3;
            if(i.first != regions[0] && key_compare(i.first,regions[0]) && regions.size() == 1){
                std::set_intersection(i.first.begin(), i.first.end(), regions[0].begin(), regions[0].end(),
                                      std::inserter(r1, r1.begin()));
                std::set_difference(regions[0].begin(), regions[0].end(),r1.begin(), r1.end(),
                                    std::inserter(r2, r2.begin()));
                std::set_difference(i.first.begin(), i.first.end(),r1.begin(), r1.end(),
                                    std::inserter(r3, r3.begin()));
                regions.emplace_back(r1);
                regions.emplace_back(r2);
                regions.emplace_back(r3);
                regions.erase(regions.begin());
            }
            else if(key_compare(i.first,regions[0]) && regions.size() > 1){
                std::set_difference(i.first.begin(), i.first.end(),regions[0].begin(), regions[0].end(),
                                    std::inserter(r3, r3.begin()));
                regions.emplace_back(r3);
            }
        }
        std::cout << "single link finished" << std::endl;
    }

    std::vector<mst_info> TLC;
    std::vector<std::set<int>> regions;
    std::map<std::set<int>, double> distance_map;
    std::map<global::pcd_vx_descriptor, global::pcd_vx_descriptor> root_m;
    std::map<global::pcd_vx_descriptor, int> root_translator;
    global::graph_v g_v;
};
