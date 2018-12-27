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
        single_linking(num_regions);
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

    void single_linking(int num_regions){
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

        for(int i=0; i<nrows; i++) printf("Region %2d: cluster %2d\n", i, clusterid[i]);

        for(int i=0; i<num_regions; i++)
            regions.emplace_back(std::set<int>());
        for(int i=0; i<nrows; i++)
            regions[clusterid[i]].insert(i);
    }

    std::vector<mst_info> TLC;
    std::vector<std::set<int>> regions;
    std::map<global::pcd_vx_descriptor, global::pcd_vx_descriptor> root_m;
    std::map<global::pcd_vx_descriptor, int> root_translator;
    global::graph_v g_v;
};
