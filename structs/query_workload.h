#ifndef EHL_QUERY_WORKLOAD_H
#define EHL_QUERY_WORKLOAD_H
#include <cmath>
#include <map>
#include <set>
#include "coverage_ordering_path.h"
#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/geometries/polygon.hpp"
#include "grid_label.h"
#include <searchinstance.h>
#include <geometry.h>
#include <ebhl.h>



namespace polyanya{
    typedef EBHL* EBHLPtr;
    const int BITS_PER_ELEMENT = sizeof(unsigned long long) * 8;


    class Workload {

    public:
        struct grid{
            polyanya::Grid_label gridLabel;
            unsigned id;
            unsigned workload;
            set<int> merged_grids;
            int section;
        };

        vector<int> mapper;
        vector<grid> grids;


        std::vector<pair<vector<Point>, vector<Point>>> visibility_area;

        Workload(int mapHeight, int mapWidth, float size, int memory, EBHLPtr eb)
                : map_height(mapHeight), map_width(mapWidth), grid_size(size), memory_limit(memory),ebhl(eb) {
            grids.clear();
            mapper.clear();

            num_of_rows = ceil((double) map_height / grid_size) ;
            num_of_columns = ceil((double) map_width / grid_size);

            grids.resize((num_of_rows)*(num_of_columns));
            for_each(grids.begin(), grids.end(), [](grid& g) {
                g.workload = 1;
            });

            traversable.resize(num_of_rows*num_of_columns);
            fill(traversable.begin(),traversable.end(),true);

            grid_label_begin.clear();
            hub_label_begin.clear();
            convex_label_set.clear();
            hub_lower_bound.clear();

        }

        Workload(int mapHeight, int mapWidth, float size, unsigned memory)
                :Workload(mapHeight, mapWidth, size, memory, nullptr) {
            num_of_rows = ceil((double) map_height / grid_size);
            num_of_columns = ceil((double) map_width / grid_size);
            grid_label_begin.clear();
            hub_label_begin.clear();
            convex_label_set.clear();
            hub_lower_bound.clear();
            grid_merge_status.clear();
            traversable.empty();
        }

        ~Workload(){
            grids.clear();
            mapper.clear();
            grid_label_begin.clear();
            hub_label_begin.clear();
            convex_label_set.clear();
            hub_lower_bound.clear();
            grid_merge_status.clear();
            traversable.empty();
        }

        void init(){
            initialise_mapper();
            for(int i = 0; i < ebhl->grid_labels.size(); ++i){
                grids[i].gridLabel = ebhl->grid_labels[i];
                grids[i].id = i;
                grids[i].merged_grids.clear();
                ebhl->grid_labels[i].hub_labels.clear();

                if(traversable[i]){
                    grid_queue.push(&grids[i]);
                }
                else{
                    grids[i].gridLabel.hub_labels.clear();
                }

                //enter mapper according to grid coordinate
                int mapper_id = floor(ebhl->grid_labels[i].a_p.x/grid_size) * num_of_rows + floor(ebhl->grid_labels[i].a_p.y/grid_size);
                mapper[mapper_id] = i;

            }
        }

        void init_parallel(){
            initialise_mapper();
            for(int i = 0; i < ebhl->grid_labels.size(); ++i){
                grids[i].gridLabel = ebhl->grid_labels[i];
                grids[i].id = i;
                grids[i].merged_grids.clear();
                ebhl->grid_labels[i].hub_labels.clear();

                //enter mapper according to grid coordinate
                int mapper_id = floor(ebhl->grid_labels[i].a_p.x/grid_size) * num_of_rows + floor(ebhl->grid_labels[i].a_p.y/grid_size);
                mapper[mapper_id] = i;


            }
        }

        void initialise_mapper(){
            //int num = max(num_of_rows,num_of_columns);
            grid_merge_status.resize(num_of_rows*num_of_columns);
            fill(grid_merge_status.begin(),grid_merge_status.end(),false);
            mapper.resize(num_of_rows*num_of_columns);
            cout << "Grids size: " << mapper.size() << endl;
            int counter = 0;
            for(int i = 0; i < num_of_rows; ++i){
                for(int j = 0; j < num_of_columns;++j){
                    mapper[counter] = -1;
                    ++counter;
                }
            }
        }

        std::vector<Point> parseLine(const std::string& line) {
            std::stringstream ss(line);
            std::string cell;
            int column = 0;
            std::vector<Point> points;
            Point point;
            bool isX = true; // Flag to track whether we are reading an x or y coordinate

            while (std::getline(ss, cell, ' ')) {
                // Read columns 4, 5 for the first point and 6, 7 for the second point
                if (column == 4 || column == 6) {
                    point.x = std::stod(cell); // Read x coordinate
                    isX = false; // Next coordinate will be y
                } else if (column == 5 || column == 7) {
                    point.y = std::stod(cell); // Read y coordinate
                    points.push_back(point); // Add the complete point to the vector
                    isX = true; // Next coordinate will be x
                }
                column++;
            }
            return points;
        }

        void load_query(const std::string& filename, int &highest_workload){
            std::ifstream file(filename);
            std::string line;
            int row = 0;

            if (!file.is_open()) {
                std::cerr << "Error: Could not open file." << std::endl;
                return;
            }

            //skip first line
            std::getline(file, line);

            while (std::getline(file, line) && row < 500000) {
                vector<Point> points = parseLine(line);
                int first_index = get_grid_index(points[0]);
                int second_index = get_grid_index(points[1]);

                grids[first_index].workload++, grids[second_index].workload++;
                if(grids[first_index].workload > highest_workload or grids[second_index].workload > highest_workload){
                    highest_workload = max(grids[first_index].workload,grids[second_index].workload);
                }

                row++;
            }
            cout << "Highest workload: " << highest_workload << " total entries read: " << row << endl;

            file.close();



        }

        int  get_grid_index(Point p) {
            return floor(p.x/grid_size) * num_of_rows + floor(p.y/grid_size);
        }

        /*
         * Find the adjacent neighbors of a given grid
         */
        set<int> find_neighbors(int grid_id) {
            set<int> neighbors;

            // Get the current grid mapper index
            int cur_grid_mapper_index = get_grid_index(grids[grid_id].gridLabel.a_p);
            vector<int> grid_to_check;
            grid_to_check.push_back(cur_grid_mapper_index);

            // Add merged grids to the list of grids to check
            if (!grids[grid_id].merged_grids.empty()) {
                for (auto id : grids[grid_id].merged_grids) {
                    grid_to_check.push_back(id);
                }
            }

            // Find neighbors for each grid in the list
            for (auto grid : grid_to_check) {

                // Check the left neighbor

                int left_grid = grid - 1;
                if(left_grid > 0){
                    int mapped_left = mapper[left_grid];
                    if(traversable[mapped_left]){
                        if (!grid_merge_status[mapped_left] && mapped_left != cur_grid_mapper_index) {
                            neighbors.insert(mapped_left);
                        }
                    }
                }
                // Check the top neighbor

                int top_grid = grid - num_of_rows;
                if(top_grid > 0){
                    int mapped_id = mapper[top_grid];

                    if(traversable[mapped_id]){
                        if (!grid_merge_status[mapped_id] && mapped_id != cur_grid_mapper_index) {
                            neighbors.insert(mapped_id);
                        }
                    }
                }


                // Check the bottom neighbor

                int bottom_grid = grid + num_of_rows;
                if(bottom_grid < grids.size()){
                    int mapped_id = mapper[bottom_grid];
                    if(traversable[mapped_id]){
                        if (!grid_merge_status[mapped_id] && mapped_id != cur_grid_mapper_index) {
                            neighbors.insert(mapped_id);
                        }
                    }
                }

                // Check the right neighbor

                int right_grid = grid + 1;
                if(right_grid < grids.size()){
                    int mapped_id = mapper[right_grid];
                    if(traversable[mapped_id]){
                        if (!grid_merge_status[mapped_id] && mapped_id != cur_grid_mapper_index) {
                            neighbors.insert(mapped_id);
                        }
                    }
                }


            }

            return neighbors;
        }

        set<int> find_neighbors_parallel(int grid_id) {
            set<int> neighbors;

            // Get the current grid mapper index
            int cur_grid_mapper_index = get_grid_index(grids[grid_id].gridLabel.a_p);
            vector<int> grid_to_check;
            grid_to_check.push_back(cur_grid_mapper_index);

            // Add merged grids to the list of grids to check
            if (!grids[grid_id].merged_grids.empty()) {
                for (auto id : grids[grid_id].merged_grids) {
                    grid_to_check.push_back(id);
                }
            }

            // Find neighbors for each grid in the list
            for (auto grid : grid_to_check) {
                int row = grid / num_of_columns;
                int col = grid % num_of_columns;

                // Check the left neighbor
                if (col > 0) {
                    int left_grid = grid - 1;
                    if (!grid_merge_status[left_grid] && left_grid != cur_grid_mapper_index && grids[left_grid].section == grids[grid_id].section) {
                        neighbors.insert(left_grid);
                    }
                }

                // Check the top neighbor
                if (row > 0) {
                    int top_grid = grid - num_of_columns;
                    if (!grid_merge_status[top_grid] && top_grid != cur_grid_mapper_index && grids[top_grid].section == grids[grid_id].section) {
                        neighbors.insert(top_grid);
                    }
                }

                // Check the bottom neighbor
                if (row < num_of_rows - 1) {
                    int bottom_grid = grid + num_of_columns;
                    if (!grid_merge_status[bottom_grid] && bottom_grid != cur_grid_mapper_index && grids[bottom_grid].section == grids[grid_id].section) {
                        neighbors.insert(bottom_grid);
                    }
                }

                // Check the right neighbor
                if (col < num_of_columns - 1) {
                    int right_grid = grid + 1;
                    if (!grid_merge_status[right_grid] && right_grid != cur_grid_mapper_index && grids[right_grid].section == grids[grid_id].section) {
                        neighbors.insert(right_grid);
                    }
                }
            }

            return neighbors;
        }

        unsigned long long get_ehl_memory(){
            unsigned long long total_memory = 0;
            for(auto gl: ebhl->grid_labels){
                for(auto label: gl.hub_labels){
                    total_memory += 4;
                    total_memory += label.convex_labels.size() * 17;
                }
            }
            return total_memory;
        }

        unsigned long long get_hub_nodes(){
            unsigned long long total_memory = 0;
            for(auto gl : grids){
                for(auto label : gl.gridLabel.hub_labels){
                    total_memory += 8;
                }
            }
            return total_memory;
        }
        unsigned long long get_total_memory(){
            unsigned long long total_memory = 0;
            for(auto gl : grids){
                for(auto label : gl.gridLabel.hub_labels){
                    total_memory += 4;
                    total_memory += label.convex_labels.size() * 17;
                }
            }
            return total_memory;
        }

        unsigned long get_memory(int grid_id){
            unsigned long memory = 0;
            for(auto label: grids[grid_id].gridLabel.hub_labels){
                memory += 4;
                memory += label.convex_labels.size() * 17;
            }
            return memory;
        }

        /*
         * Computation of similarity score of adjacent cells identified
         */
        vector<pair<int,double>> compute_neighbor_jaccard_similarity(set<int>& neighbors, int grid_id){
            vector<pair<int,double>> similarity;

            set<int> grid_hubs;
            for(const auto &hub: grids[grid_id].gridLabel.hub_labels){
                grid_hubs.insert(hub.hub_id);
            }

            for(int i : neighbors){
                set<int> neighbor_hubs;
                int neighbor_id = mapper[i];
                for(const auto &hub: grids[neighbor_id].gridLabel.hub_labels){
                    neighbor_hubs.insert(hub.hub_id);
                }

                //calcualte intersection
                vector<int> intersection;
                set_intersection(grid_hubs.begin(), grid_hubs.end(), neighbor_hubs.begin(), neighbor_hubs.end(),
                                 std::back_inserter(intersection));

                // Calculate union
                std::vector<int> union_set;
                std::set_union(grid_hubs.begin(), grid_hubs.end(), neighbor_hubs.begin(), neighbor_hubs.end(),
                               std::back_inserter(union_set));

                if (union_set.empty()) similarity.push_back(make_pair(i,0));  // Avoid division by zero
                else{

                    //add in workload factor
                    double workload = grids[neighbor_id].workload;
                    if(grids[neighbor_id].workload == 0){
                        workload = 0.5;
                    }
                    double workload_factor = (1.0/workload) * 0.2;
                    double sim_factor = (static_cast<double>(intersection.size()) / static_cast<double>(union_set.size())) * 0.8;


                    double score = sim_factor + workload_factor;
                    //cout << "Score: " << score << " workload: " << workload_factor << endl;


                    similarity.push_back(make_pair(i,score));
                }

            }

            sort(similarity.begin(), similarity.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second > b.second;
            });

            return similarity;
        }


        /*
         * Process of merging one grid into another
         */
        void merge_grids(int grid_id, int neighbor_id) {
            grids[grid_id].workload += grids[neighbor_id].workload;

            assert(!grids[grid_id].gridLabel.hub_labels.empty());
            assert(!grids[neighbor_id].gridLabel.hub_labels.empty());

            // Use a map to merge hub labels from both grids and a set to track updated hubs
            std::map<int, polyanya::Hub_label> label_mapper;
            std::set<int> updated_hubs;

            // Add labels from the first grid
            for (auto &hub : grids[grid_id].gridLabel.hub_labels) {
                label_mapper[hub.hub_id] = hub;
            }

            // Merge labels from the second grid
            for (auto &hub : grids[neighbor_id].gridLabel.hub_labels) {
                if (label_mapper.find(hub.hub_id) == label_mapper.end()) {
                    // If hub_id not found, add the whole hub label with false visibility
                    polyanya::Hub_label new_hub_label = hub;
                    for (auto &convex : new_hub_label.convex_labels) {
                        convex.visibility = false;
                    }
                    label_mapper[hub.hub_id] = new_hub_label;
                } else {
                    // Hub exists, mark it as updated
                    updated_hubs.insert(hub.hub_id);

                    label_mapper[hub.hub_id].min_lower_bound = std::min(label_mapper[hub.hub_id].min_lower_bound, hub.min_lower_bound);

                    // Create a set to track updated convex labels
                    std::set<int> updated_convex;

                    // Merge convex labels
                    auto &existing_labels = label_mapper[hub.hub_id].convex_labels;
                    for (const auto &convex : hub.convex_labels) {
                        auto it = std::find_if(existing_labels.begin(), existing_labels.end(),
                                               [&convex](const Convex_vertices_label &item) {
                                                   return item.convex_vertex == convex.convex_vertex;
                                               });

                        if (it != existing_labels.end()) {
                            // Existing label found
                            // Update visibility if necessary
                            if (!it->visibility or !convex.visibility) {
                                it->visibility = false;
                            }
                            updated_convex.insert(it->convex_vertex);  // Mark this convex as updated
                        } else {
                            // Add new label with visibility set to false, as it does not exist
                            Convex_vertices_label new_convex = convex;
                            new_convex.visibility = false;  // Set visibility to false for new convex label
                            existing_labels.push_back(new_convex);
                        }
                    }

                    // Set visibility to false for any convex label not updated from the second grid
                    for (auto &convex : existing_labels) {
                        if (updated_convex.find(convex.convex_vertex) == updated_convex.end()) {
                            convex.visibility = false;
                        }
                    }
                }
            }

            // Set visibility to false for any hub not updated from the second grid
            for (auto &pair : label_mapper) {
                if (updated_hubs.find(pair.first) == updated_hubs.end()) {
                    // This hub was not updated, set all convex label visibilities to false
                    for (auto &convex : pair.second.convex_labels) {
                        convex.visibility = false;
                    }
                }
            }

            // Clear the old hub labels
            grids[grid_id].gridLabel.hub_labels.clear();

            // Convert the map back to the vector
            for (const auto &pair : label_mapper) {
                grids[grid_id].gridLabel.hub_labels.push_back(pair.second);
            }

            // Sort the hub labels by hub_id
            std::sort(grids[grid_id].gridLabel.hub_labels.begin(), grids[grid_id].gridLabel.hub_labels.end(),
                      [](const polyanya::Hub_label &l, const polyanya::Hub_label &r) {
                          return l.hub_id < r.hub_id;
                      });
        }

        /*
         * Updates the mapper of merged grids with the new region it is merged into
         */
        void update_mapper(int neighbor_grid_id, int new_id){
            //iterates over all the grids in the region to be merged
            for(auto grid: grids[neighbor_grid_id].merged_grids){
                //inserts the grids into the region being merged and updates the mapper
                grids[new_id].merged_grids.insert(grid);
                mapper[grid] = new_id;
            }

            grids[neighbor_grid_id].merged_grids.clear();
        }

        /*
         * Empties and rebuild the queue every 1000 iteration
         */
        void clean_up_queue(){
            pq empty_queue;
            swap(grid_queue,empty_queue);

            for(int i = 0; i < grids.size(); ++i){
                if(traversable[i]){
                    if(!grid_merge_status[i]){
                        grid_queue.push(&grids[i]);
                    }
                }
            }

        }

        void build_queue_parallel(int num_threads){
            for(int i = 0; i < grid_queue_parallel.size(); ++i){
                pq empty_queue;
                swap(grid_queue_parallel[i],empty_queue);
            }

            int sections_x = sqrt(num_threads);
            int sections_y = sqrt(num_threads);
            int total_sections = sections_x * sections_y;

            float section_width = static_cast<float>(map_width) / sections_y;
            float section_height = static_cast<float>(map_height) / sections_x;

            grid_queue_parallel.resize(9);


            for(int i = 0; i < grids.size(); ++i){

                // Determine which section the grid cell belongs to
                int grid_x = static_cast<int>(grids[i].gridLabel.a_p.x);
                int grid_y = static_cast<int>(grids[i].gridLabel.a_p.y);
                int section_x = grid_x / section_width;
                int section_y = grid_y / section_height;
                int section_id = section_y * sections_x + section_x;

                // Ensure section_id is within bounds
                if (section_id >= total_sections) {
                    std::cerr << "Section ID out of bounds: " << section_id << std::endl;
                }

                grids[i].section = section_id;
                grid_queue_parallel[section_id].push(&grids[i]);

            }
        }

        /*
         * Compression function for EHL*
         * Detailed description is in Algorithm 1 in the paper
         */
        void merge(){
            unsigned long long total_memory = get_total_memory();

            memory_limit = ceil(((double)memory_limit/100) * total_memory);
            cout << "Total memory: " << total_memory << " ,Memory limit: " <<  memory_limit << endl;

            int operation_count = 0;
            while(total_memory > memory_limit){
                operation_count++;
                if(grid_queue.empty()){
                    break;
                }
                auto cur_grid = grid_queue.top();
                grid_queue.pop();

                //check if it's been merged before
                if(grid_merge_status[cur_grid->id]){
                    continue;
                }
                //find the four neighbors of current grid
                set<int> grid_neighbors = find_neighbors(cur_grid->id);

                //compute jaccard similarity score for the four neighbors
                vector<pair<int,double>> similarity = compute_neighbor_jaccard_similarity(grid_neighbors,cur_grid->id);


                if(similarity.empty() ){
                    continue;
                }
                //merge the current grid with the neighbor having highest jaccard score
                int neighbor_grid_id = mapper[similarity[0].first];

                //update memory
                total_memory -= get_memory(cur_grid->id);
                total_memory -= get_memory(neighbor_grid_id);

                //merges the two grids
                merge_grids(cur_grid->id, neighbor_grid_id);

                //store merged cell as mapper id
                grids[cur_grid->id].merged_grids.insert(similarity[0].first);

                //as neighbor is merged, status is changed
                grid_merge_status[neighbor_grid_id] = true;
                grids[neighbor_grid_id].gridLabel.hub_labels.clear();

                //update mapper of merged neighbor grid, if merged grid is a region, all the grids inside the region has the mapper updated
                mapper[similarity[0].first] = cur_grid->id;
                update_mapper(neighbor_grid_id,cur_grid->id);

                grid_queue.push(&grids[cur_grid->id]);


                //update memory
                total_memory += get_memory(cur_grid->id);

                //cleans up the queue after every 1000 iterations to clean out all grids that has already been merged to speed up performance
                if(operation_count == 1000){
                    clean_up_queue();
                    operation_count = 0;
                }

            }

            cout << total_memory << " " << get_total_memory() << endl;

        }

        double get_avg_grid(){
            int merged_grids = 0;
            int grids_size = 0;
            for(int i = 0; i < grid_merge_status.size(); ++i){
                if(traversable[i]){
                    if(!grid_merge_status[i]){
                        merged_grids++;
                        grids_size += grids[i].merged_grids.size()+1;
                    }
                }

            }
            return grids_size/(double)merged_grids;
        }


        void save_adjacent_workload(const char *save_filename){
            cout << "Saving workload..........." << endl;
            unsigned g_label_start = 0;
            unsigned c_label_start = 0;
            for(int i = 0; i < grids.size(); ++i){
                grid_label_begin.push_back(g_label_start);
                for(int j = 0; j < grids[i].gridLabel.hub_labels.size(); ++j){
                    g_label_start++;
                    g_label_start++;
                    hub_label_begin.push_back(grids[i].gridLabel.hub_labels[j].hub_id);
                    hub_label_begin.push_back(c_label_start);
                    hub_lower_bound.push_back(grids[i].gridLabel.hub_labels[j].min_lower_bound);
                    for(int k = 0; k < grids[i].gridLabel.hub_labels[j].convex_labels.size(); ++k){
                        convex_label_set.push_back(grids[i].gridLabel.hub_labels[j].convex_labels[k]);
                        c_label_start++;
                    }
                }
            }
            hub_label_begin.push_back(0);
            hub_label_begin.push_back(c_label_start);
            grid_label_begin.push_back(g_label_start);
            hub_label_begin.shrink_to_fit();
            grid_label_begin.shrink_to_fit();
            convex_label_set.shrink_to_fit();
            hub_lower_bound.shrink_to_fit();
            string filename = save_filename;
            ofstream ofs(filename, ios::binary | ios::out);
            int size = grids.size();
            ofs.write((const char *) &size, sizeof(size));
            for (int i = 0; i < size; ++i) {
                int hub_size = grids[i].gridLabel.hub_labels.size();
                ofs.write((const char *) &hub_size, sizeof(hub_size));
                ofs.write((const char *) &grids[i].gridLabel.a_p, sizeof(grids[i].gridLabel.a_p));
                ofs.write((const char *) &grids[i].gridLabel.b_p, sizeof(grids[i].gridLabel.b_p));
                ofs.write((const char *) &grids[i].gridLabel.c_p, sizeof(grids[i].gridLabel.c_p));
                ofs.write((const char *) &grids[i].gridLabel.d_p, sizeof(grids[i].gridLabel.d_p));
                ofs.write((const char *) &grids[i].id, sizeof(grids[i].id));
            }

            cout << "save hub label begin" << endl;
            int hub_label_size = hub_label_begin.size();
            ofs.write((const char *) &hub_label_size, sizeof(hub_label_size));
            for (int i = 0; i < hub_label_size; ++i) {
                ofs.write((const char *) &hub_label_begin[i], sizeof(hub_label_begin[i]));
            }
            int grid_label_size = grid_label_begin.size();
            ofs.write((const char *) &grid_label_size, sizeof(grid_label_size));
            for (int i = 0; i < grid_label_size; ++i) {
                ofs.write((const char *) &grid_label_begin[i], sizeof(grid_label_begin[i]));
            }
            int convex_vertex_size = convex_label_set.size();
            ofs.write((const char *) &convex_vertex_size, sizeof(convex_vertex_size));
            for (int i = 0; i < convex_vertex_size; ++i) {
                ofs.write((const char *) &convex_label_set[i].convex_vertex, sizeof(convex_label_set[i].convex_vertex));
                ofs.write((const char *) &convex_label_set[i].distance, sizeof(convex_label_set[i].distance));
                ofs.write((const char *) &convex_label_set[i].predecessor, sizeof(convex_label_set[i].predecessor));
                ofs.write((const char *) &convex_label_set[i].visibility, sizeof(convex_label_set[i].visibility));
            }

            int hub_lower_size = hub_lower_bound.size();
            ofs.write((const char *) &hub_lower_size, sizeof(hub_lower_size));
            for(int i = 0; i < hub_lower_size; ++i){
                ofs.write((const char *) &hub_lower_bound[i], sizeof(hub_lower_bound[i]));
            }
            ofs.close();
        }

        void load_adjacent_workload(const char *load_filename) {
            cout << "loading workload grids..." << endl;
            ifstream ifs(load_filename, ios::binary | ios::in);
            if(!ifs.good()){
                cout << "cannot read file" << endl;
            }
            int isize = 0;
            ifs.read((char *) &isize, sizeof(isize));
            grids.resize(isize);
            cout << "grids: " << grids.size() << endl;
            for (int i = 0; i < isize; ++i) {
                int hub_size = 0;
                ifs.read((char *) &hub_size, sizeof(hub_size));
                grid &node = grids[i];
                node.gridLabel.hub_labels.resize(hub_size);
                ifs.read((char *) &node.gridLabel.a_p, sizeof(node.gridLabel.a_p));
                ifs.read((char *) &node.gridLabel.b_p, sizeof(node.gridLabel.b_p));
                ifs.read((char *) &node.gridLabel.c_p, sizeof(node.gridLabel.c_p));
                ifs.read((char *) &node.gridLabel.d_p, sizeof(node.gridLabel.d_p));
                ifs.read((char *) &node.id, sizeof(node.id));
            }
            int hub_label_size;
            ifs.read((char *) &hub_label_size, sizeof(hub_label_size));
            hub_label_begin.resize(hub_label_size);
            for(int i = 0; i < hub_label_size; ++i){
                unsigned &h = hub_label_begin[i];
                ifs.read((char *) &h, sizeof(h));
            }
            int grid_label_size;
            ifs.read((char *) &grid_label_size, sizeof(grid_label_size));
            grid_label_begin.resize(grid_label_size);
            for (int i = 0; i < grid_label_size; ++i) {
                unsigned &g = grid_label_begin[i];
                ifs.read((char *) &g, sizeof(g));
            }

            int convex_vertex_size;
            ifs.read((char *) &convex_vertex_size, sizeof(convex_vertex_size));
            convex_label_set.resize(convex_vertex_size);
            for (int i = 0; i < convex_vertex_size; ++i) {
                Convex_vertices_label &label = convex_label_set[i];
                ifs.read((char *) &label.convex_vertex, sizeof(label.convex_vertex));
                ifs.read((char *) &label.distance, sizeof(label.distance));
                ifs.read((char *) &label.predecessor, sizeof(label.predecessor));
                ifs.read((char *) &label.visibility, sizeof(label.visibility));
            }
            int hub_lower_size;
            ifs.read((char *) &hub_lower_size, sizeof(hub_lower_size));
            hub_lower_bound.resize(hub_lower_size);
            for(int i = 0; i < hub_lower_size; ++i){
                double &lower_bound = hub_lower_bound[i];
                ifs.read((char *) &lower_bound, sizeof(lower_bound));
            }
            hub_label_begin.shrink_to_fit();
            grid_label_begin.shrink_to_fit();
            convex_label_set.shrink_to_fit();
            hub_lower_bound.shrink_to_fit();
            ifs.close();
        }
        void load_non_taut_triangles(const char *filename) {
            ifstream in(filename);
            vector<vector<double>> lines;
            int size = 0;
            if (in) {
                string line;
                getline(in, line);
                size = stoi(line);
                while (getline(in, line)) {
                    stringstream sep(line);
                    string field;
                    lines.push_back(vector<double>());
                    while (getline(sep, field, ' ')) {
                        lines.back().push_back(stod(field));
                    }
                }
            }
            in.close();
            visibility_area.resize(size);
            for (auto line : lines) {
                int vertex = line[0];
                int poly_id = vertex / 2;
                if (vertex % 2 == 0) {
                    visibility_area[poly_id].first.push_back(Point{line[1], line[2]});
                } else {
                    visibility_area[poly_id].second.push_back(Point{line[1], line[2]});
                }
            }
            for(auto &va :  visibility_area){
                va.first.erase(remove(va.first.begin(),va.first.end(), va.first.front())-1,va.first.end());
                va.first.shrink_to_fit();
                va.second.erase(remove(va.second.begin(),va.second.end(), va.second.front())-1,va.second.end());
                va.second.shrink_to_fit();
            }
            visibility_area.shrink_to_fit();
            std::cout<<"Finish loading non-taut visible area"<<std::endl;
        }

        void load_quad_mapper(const char *filename){

            std::ifstream inputFile(filename, std::ios::binary);
            if (!inputFile) {
                std::cout << "Failed to open the file for reading." << std::endl;
            }
            size_t size;
            inputFile.read(reinterpret_cast<char*>(&size), sizeof(size));
            mapper.resize(size);
            inputFile.read(reinterpret_cast<char*>(mapper.data()), size * sizeof(int));

            inputFile.close();
        }

        const grid* get_grid(Point p){
            int mapper_id = floor(p.x/grid_size) * num_of_rows + floor(p.y/grid_size);
            grid *cur_node = &grids[mapper[mapper_id]];
            return cur_node;
        }

        const unsigned& get_grid_labels_begin(const unsigned& grid_index){
            assert(grid_index <= grid_label_begin.size());
            return grid_label_begin[grid_index];
        }

        const unsigned& get_grid_labels_end(const unsigned& grid_index){
            assert(grid_index +1 <= grid_label_begin.size());
            return grid_label_begin[grid_index+1];
        }

        const unsigned& get_hub_id(const unsigned& hub_label_index){
            return hub_label_begin[hub_label_index];
        }

        const double& get_lower_bound(const unsigned& hub_label_index){
//            assert(hub_label_index*2+1 <= hub_label_begin.size());
            return hub_lower_bound[hub_label_index/2];
        }

        const unsigned& get_convex_begin(const unsigned& hub_label_index){
//            assert(hub_label_index*2+1 <= hub_label_begin.size());
            return hub_label_begin[hub_label_index + 1];
        }

        const unsigned& get_convex_end(const unsigned& hub_label_index){
            //assert((hub_label_index+1)*2+1 <= hub_label_begin.size());
            return hub_label_begin[hub_label_index+ 3];
        }

        const Convex_vertices_label& get_convex_label(const unsigned& convex_label_index){
            return convex_label_set[convex_label_index];
        }

        static unsigned short get_orientation_value(
                const Point& a, const Point& b, const Point& c
        )
        {
            const double cross = (b - a) * (c - b);
            if (std::abs(cross) < EPSILON)
            {
//                return Orientation::COLLINEAR;
                return 1;
            } else if (cross > 0)
            {
//                return Orientation::CCW;
                return 2 ;
            }
            else
            {
//                return Orientation::CW;
                return 0;
            }
        }


        bool check_visiblity(const Point& check_point, const Point& obstacle_middle_point,const Point& turning_point, unsigned int convex_vertex_id){
            assert(convex_vertex_id < visibility_area.size());
            const pair<vector<Point>,vector<Point> >& va = visibility_area[convex_vertex_id];
            assert(va.first.size() >= 2);
            const Orientation& o1 = get_orientation(obstacle_middle_point, turning_point,va.first[va.first.size()/2] );
            const Orientation& o2 = get_orientation(obstacle_middle_point, turning_point, check_point);
            if(o1 == o2){
                const Orientation& check_ori1 = get_orientation( turning_point, va.first.front(),check_point);
                const Orientation& check_ori2 = get_orientation( turning_point, va.first.back(), check_point);
                if(check_ori1 ==  Orientation::COLLINEAR){
                    return onSegment(turning_point, check_point,va.first.front());

                }else if( check_ori2 ==  Orientation::COLLINEAR){
                    return onSegment(turning_point, check_point,va.first.back());

                }else if( check_ori1 == Orientation::CW && check_ori2 == Orientation::CCW ){
                    const auto&  it = std::lower_bound(va.first.begin(), va.first.end(), 1,
                                                       [&turning_point, &check_point](const Point& lvalue, const unsigned short& value){
                                                           return get_orientation_value(turning_point, lvalue,check_point) < value;
//                                        double angle = get_relative_angle(obstacle_middle_point,turning_point,lvalue) + EPSILON;
//                                         return angle < value;
                                                       });

                    if(is_collinear(turning_point,*it,check_point)){
                        const Orientation temp_ori1 = get_orientation(*(it-1), check_point, *it);
                        if(temp_ori1 != Orientation :: CW){
                            if(temp_ori1 == Orientation::COLLINEAR){
                                return onSegment(*it, check_point, *(it-1));
                            }
                            else{
                                return true;
                            }

//                        if(get_orientation(  *(it-1), check_point, *it) != Orientation ::CW){
//                            return true;
                        }else{
                            const Orientation& ori2 = get_orientation( *it, check_point,*(it+1));
                            if(ori2 == Orientation::COLLINEAR){
                                return onSegment(*it, check_point,*(it+1));
                            }else{
                                return  ori2 != Orientation ::CW;
                            }
                        }
                    }else{
                        return get_orientation( *(it-1), check_point, *it) != Orientation ::CW;
                    }

                }else{
                    return false;
                }


            }else{
                const Orientation& check_ori1 = get_orientation( turning_point, va.second.front(),check_point);
                const Orientation& check_ori2 = get_orientation( turning_point, va.second.back(), check_point);
                if(check_ori1 ==  Orientation::COLLINEAR){
                    return onSegment(turning_point, check_point,va.second.front());

                }else if( check_ori2 ==  Orientation::COLLINEAR){
                    return onSegment(turning_point, check_point,va.second.back());

                }else if( check_ori1 == Orientation::CW && check_ori2 == Orientation::CCW ){
//                    const auto&  it = std::lower_bound(va.second.begin(), va.second.end(), check_point_angle,
//                                                       [&turning_point, &obstacle_middle_point](const Point& lvalue, double value){
//                                                           double angle = get_relative_angle(obstacle_middle_point,turning_point,lvalue) + EPSILON;
//                                                           return angle < value;
//                                                       });
                    const auto&  it = std::lower_bound(va.second.begin(), va.second.end(), 1,
                                                       [&turning_point, &check_point](const Point& lvalue, const unsigned short& value){
                                                           return get_orientation_value(turning_point, lvalue,check_point) < value;
                                                       });

                    if(is_collinear(turning_point,*it,check_point)){
                        const Orientation temp_ori = get_orientation(*(it-1), check_point, *it);
                        if(temp_ori != Orientation :: CW){
                            if(temp_ori == Orientation::COLLINEAR){
                                return onSegment(*it, check_point, *(it-1));
                            }
                            else{
                                return true;
                            }

                            //if(get_orientation(  *(it-1), check_point, *it) != Orientation ::CW){
                            //return true;
                        }else{
                            const Orientation& ori2 = get_orientation( *it, check_point,*(it+1));
                            if(ori2 == Orientation::COLLINEAR){
                                return onSegment(*it, check_point,*(it+1));
                            }else{
                                return  ori2 != Orientation ::CW;
                            }
                        }
                    }else{
                        return get_orientation( *(it-1), check_point, *it) != Orientation ::CW;
                    }

                }else{
                    return false;
                }
            }
        }


    //private:

        typedef grid* gridPtr;

        int map_height,map_width,num_of_rows,num_of_columns;
        float grid_size;
        unsigned long long memory_limit;
        const EBHLPtr ebhl;

        std::vector<unsigned> grid_label_begin;
        std::vector<unsigned> hub_label_begin;
        std::vector<double> hub_lower_bound;
        std::vector<Convex_vertices_label> convex_label_set;

        struct CompareWorkload {
            bool operator()(const gridPtr& a, const gridPtr& b) {
                return a->workload > b->workload;  // This ensures the grid with the lowest workload is at the top
            }
        };
        typedef priority_queue<gridPtr,vector<gridPtr>,CompareWorkload> pq;
        pq grid_queue;
        vector<bool> grid_merge_status;
        vector<bool> traversable;

        vector<pq> grid_queue_parallel;
    };
}



#endif //EHL_QUERY_WORKLOAD_H
