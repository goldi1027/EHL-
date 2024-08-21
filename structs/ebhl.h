/*
 * Implementation of EHL structure
 * Includes initialisation of different versions of EHL, checking visibility and I/O
 */

#ifndef LEARNED_HEURISTIC_EBHL_H
#define LEARNED_HEURISTIC_EBHL_H
#include <grid_label.h>
#include <ostream>
#include <iostream>
#include <fstream>
#include "searchinstance.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <polygon.h>
namespace bg = boost::geometry;
namespace polyanya{
    class EBHL {
        typedef boost::geometry::model::d2::point_xy<double> point_xy;
        typedef boost::geometry::model::polygon<point_xy> polygon;
        typedef boost::geometry::strategy::within::franklin<point_xy, point_xy, void> fran;


    public:
        std::vector<unsigned> grid_label_begin;
        std::vector<unsigned> hub_label_begin;
        std::vector<double> hub_lower_bound;

        std::vector<Grid_label> grid_labels;
        std::vector<pair<vector<Point>, vector<Point>>> visibility_area;
        std::vector<Convex_vertices_label> convex_label_set;

        int grid_size;
        int map_height;
        int map_width;
        int num_of_columns;
        int num_of_rows;
        int current_grid_size;
        vector<Point> grid_corner_points;

        EBHL(int gridSize, int mapHeight, int mapWidth) {
            grid_size = gridSize;
            map_height = mapHeight;
            map_width = mapWidth;

            num_of_rows = ceil((double) map_height / grid_size);
            num_of_columns = ceil((double) map_width / grid_size);
            grid_label_begin.clear();
            hub_label_begin.clear();
            hub_lower_bound.clear();
            convex_label_set.clear();
            visibility_area.clear();
        };

        ~EBHL() {
            grid_label_begin.clear();
            hub_label_begin.clear();
            convex_label_set.clear();
            visibility_area.clear();
            hub_lower_bound.clear();
        };

        /*
         * Initialise EHL instance
         * Superimpose grids with a given grid size for the entire map
         */
        void initialize_grid_map(int grid_size, int map_width, int map_height) {
            int row = ceil((double) map_width / grid_size);
            int column = ceil((double) map_height / grid_size);
            int map_size = row * column;
            grid_labels.resize(map_size);
            int counter = 0;
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < column; ++j) {
                    grid_labels[counter].a_p.x = i * grid_size;
                    grid_labels[counter].a_p.y = j * grid_size;
                    grid_labels[counter].b_p.x = i * grid_size + grid_size;
                    grid_labels[counter].b_p.y = j * grid_size;
                    grid_labels[counter].c_p.x = i * grid_size + grid_size;;
                    grid_labels[counter].c_p.y = j * grid_size + grid_size;;
                    grid_labels[counter].d_p.x = i * grid_size;
                    grid_labels[counter].d_p.y = j * grid_size + grid_size;;
                    ++counter;
                }
            }

        }

        void save_grid_labels_parallel(const char *save_filename) {
            cout << "Saving grid labels in parallel..." << endl;
            {
                printf("Using %d threads\n", omp_get_max_threads());
#pragma omp parallel
                {
                    const int thread_count = omp_get_num_threads();
                    const int thread_id = omp_get_thread_num();
                    const int node_count = grid_labels.size();
                    int node_begin = (node_count * thread_id) / thread_count;
                    int node_end = (node_count * (thread_id + 1)) / thread_count;
                    string file_name = save_filename;
                    string filename = file_name + "_" + to_string(thread_id);
//                    cout << thread_id << endl;
                    ofstream ofs(filename, ios::binary | ios::out);
                    int size = grid_labels.size();
                    ofs.write((const char *) &size, sizeof(size));
                    for (int source_node = node_begin; source_node < node_end; ++source_node) {
                        int visible_size = grid_labels[source_node].hub_labels.size();
                        //int p_visible_size = grid_labels[v].p_hub_labels.size();

                        ofs.write((const char *) &visible_size, sizeof(visible_size));
                        //ofs.write((const char*)&p_visible_size, sizeof(p_visible_size));

                        ofs.write((const char *) &grid_labels[source_node].a_p, sizeof(grid_labels[source_node].a_p));
                        ofs.write((const char *) &grid_labels[source_node].b_p, sizeof(grid_labels[source_node].b_p));
                        ofs.write((const char *) &grid_labels[source_node].c_p, sizeof(grid_labels[source_node].c_p));
                        ofs.write((const char *) &grid_labels[source_node].d_p, sizeof(grid_labels[source_node].d_p));

                        for (int i = 0; i < grid_labels[source_node].hub_labels.size(); ++i) {
                            int temp_size = grid_labels[source_node].hub_labels[i].convex_labels.size();
                            ofs.write((const char *) &grid_labels[source_node].hub_labels[i].hub_id,
                                      sizeof(grid_labels[source_node].hub_labels[i].hub_id));
                            ofs.write((const char *) &temp_size, sizeof(temp_size));
                            for (int j = 0; j < temp_size; ++j) {
                                ofs.write(
                                        (const char *) &grid_labels[source_node].hub_labels[i].convex_labels[j].convex_vertex,
                                        sizeof(grid_labels[source_node].hub_labels[i].convex_labels[j].convex_vertex));
                                ofs.write(
                                        (const char *) &grid_labels[source_node].hub_labels[i].convex_labels[j].distance,
                                        sizeof(grid_labels[source_node].hub_labels[i].convex_labels[j].distance));
                                ofs.write(
                                        (const char *) &grid_labels[source_node].hub_labels[i].convex_labels[j].visibility,
                                        sizeof(grid_labels[source_node].hub_labels[i].convex_labels[j].visibility));
                                ofs.write(
                                        (const char *) &grid_labels[source_node].hub_labels[i].convex_labels[j].predecessor,
                                        sizeof(grid_labels[source_node].hub_labels[i].convex_labels[j].predecessor));
                            }
                            ofs.write((const char *) &grid_labels[source_node].hub_labels[i].min_lower_bound,
                                      sizeof(grid_labels[source_node].hub_labels[i].min_lower_bound));
                        }
                    }
                    ofs.close();
                }
            }
        }

        void save_adjacent_list(const char *save_filename) {
            std::vector<unsigned> grid_label_begin;
            std::vector<unsigned> hub_label_begin;
            std::vector<Convex_vertices_label> convex_label_set;
            vector<double> hub_lower_bound;
            unsigned g_label_start = 0;
            unsigned c_label_start = 0;
            for (const auto &g_label: grid_labels) {
                grid_label_begin.push_back(g_label_start);
                for (const auto &h_label: g_label.hub_labels) {
                    g_label_start++;
                    g_label_start++;
                    hub_label_begin.push_back(h_label.hub_id);
                    hub_label_begin.push_back(c_label_start);
                    hub_lower_bound.push_back(h_label.min_lower_bound);
                    for (const auto &c_label: h_label.convex_labels) {
                        convex_label_set.push_back(c_label);
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

        void load_adjacent_list(const char *load_filename) {
            cout << "Loading adjacent list " << endl;
            string file_name = load_filename;
            ifstream ifs(file_name);
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
            ifs.close();
            hub_label_begin.shrink_to_fit();
            grid_label_begin.shrink_to_fit();
            convex_label_set.shrink_to_fit();
            hub_lower_bound.shrink_to_fit();
        }

        void save_grid_labels(const char *save_filename) {
            cout << "Saving grid labels!" << endl;
            ofstream ofs(save_filename, ios::binary | ios::out);
            int size = grid_labels.size();
            ofs.write((const char *) &size, sizeof(size));
            unsigned long long total_number_of_label = 0;
            for (int v = 0; v < grid_labels.size(); ++v) {
                double ratio = (double) v / grid_labels.size() * 100;
                std::cout << "Progress: [" << v << "/" << grid_labels.size() << "] "
                          << std::setprecision(3) << ratio << "% \r";
                std::cout.flush();

                int visible_size = grid_labels[v].hub_labels.size();
                //int p_visible_size = grid_labels[v].p_hub_labels.size();
                ofs.write((const char *) &visible_size, sizeof(visible_size));
                //ofs.write((const char*)&p_visible_size, sizeof(p_visible_size));

                ofs.write((const char *) &grid_labels[v].a_p, sizeof(grid_labels[v].a_p));
                ofs.write((const char *) &grid_labels[v].b_p, sizeof(grid_labels[v].b_p));
                ofs.write((const char *) &grid_labels[v].c_p, sizeof(grid_labels[v].c_p));
                ofs.write((const char *) &grid_labels[v].d_p, sizeof(grid_labels[v].d_p));

                for (int i = 0; i < grid_labels[v].hub_labels.size(); ++i) {
                    int temp_size = grid_labels[v].hub_labels[i].convex_labels.size();
                    total_number_of_label += temp_size;
                    ofs.write((const char *) &grid_labels[v].hub_labels[i].hub_id,
                              sizeof(grid_labels[v].hub_labels[i].hub_id));
                    ofs.write((const char *) &temp_size, sizeof(temp_size));
                    for (int j = 0; j < temp_size; ++j) {
                        ofs.write((const char *) &grid_labels[v].hub_labels[i].convex_labels[j].convex_vertex,
                                  sizeof(grid_labels[v].hub_labels[i].convex_labels[j].convex_vertex));
                        ofs.write((const char *) &grid_labels[v].hub_labels[i].convex_labels[j].distance,
                                  sizeof(grid_labels[v].hub_labels[i].convex_labels[j].distance));
                        ofs.write((const char *) &grid_labels[v].hub_labels[i].convex_labels[j].visibility,
                                  sizeof(grid_labels[v].hub_labels[i].convex_labels[j].visibility));
                        ofs.write((const char *) &grid_labels[v].hub_labels[i].convex_labels[j].predecessor,
                                  sizeof(grid_labels[v].hub_labels[i].convex_labels[j].predecessor));
                    }
                    ofs.write((const char *) &grid_labels[v].hub_labels[i].min_lower_bound,
                              sizeof(grid_labels[v].hub_labels[i].min_lower_bound));
                }
            }
            std::cout << "Total number of label:" << total_number_of_label << std::endl;
            std::cout << "Average label per grid:" << (double) total_number_of_label / size << std::endl;
            ofs.close();
        }

        void load_grid_labels(const char *filename) {
            ifstream ifs(filename);
            if (!ifs) {
                perror("Failed because: ");
                exit(1);
            }

            int isize = 0;
            ifs.read((char *) &isize, sizeof(isize));
            grid_labels.resize(isize);

            for (int v = 0; v < isize; ++v) {
                Grid_label &l = grid_labels[v];
                int visible_size;
                ifs.read((char *) &visible_size, sizeof(visible_size));
                l.hub_labels.resize(visible_size);
                ifs.read((char *) &l.a_p, sizeof(l.a_p));
                ifs.read((char *) &l.b_p, sizeof(l.b_p));
                ifs.read((char *) &l.c_p, sizeof(l.c_p));
                ifs.read((char *) &l.d_p, sizeof(l.d_p));
                if (visible_size > 0) {
                    for (int i = 0; i < visible_size; ++i) {
                        ifs.read((char *) &l.hub_labels[i].hub_id, sizeof(l.hub_labels[i].hub_id));
                        int convex_size;
                        ifs.read((char *) &convex_size, sizeof(convex_size));
                        l.hub_labels[i].convex_labels.resize(convex_size);
                        for (int j = 0; j < convex_size; ++j) {
                            ifs.read((char *) &l.hub_labels[i].convex_labels[j].convex_vertex,
                                     sizeof(l.hub_labels[i].convex_labels[j].convex_vertex));
                            ifs.read((char *) &l.hub_labels[i].convex_labels[j].distance,
                                     sizeof(l.hub_labels[i].convex_labels[j].distance));
                            ifs.read((char *) &l.hub_labels[i].convex_labels[j].visibility,
                                     sizeof(l.hub_labels[i].convex_labels[j].visibility));
                            ifs.read((char *) &l.hub_labels[i].convex_labels[j].predecessor,
                                     sizeof(l.hub_labels[i].convex_labels[j].predecessor));
                        }
                        ifs.read((char *) &l.hub_labels[i].min_lower_bound,
                                 sizeof(l.hub_labels[i].min_lower_bound));
                    }
                }
            }
            ifs.close();
            std::cout << "Finish loading grid label" << std::endl;
        }


        void load_grid_labels_parallel(const char *load_filename) {
            cout << "loading grid labels in parallel..." << endl;
            //load grid label size of the first grid label file
            string file_name = load_filename;
            string first_file = file_name + "_0";
            ifstream ifs(first_file);
            int isize = 0;
            ifs.read((char *) &isize, sizeof(isize));
            grid_labels.resize(isize);
            {
                printf("Using %d threads\n", omp_get_max_threads());
#pragma omp parallel
                {
                    const int thread_count = omp_get_num_threads();
                    const int thread_id = omp_get_thread_num();
                    const int node_count = grid_labels.size();
                    int node_begin = (node_count * thread_id) / thread_count;
                    int node_end = (node_count * (thread_id + 1)) / thread_count;
                    string filename = file_name + "_" + to_string(thread_id);
//                    cout << thread_id << endl;
                    ifstream ifs(filename);
                    cout << filename << endl;
                    //because all the grid label files will write the size of grid label, we could just throw it away and do nothing with it
                    int temp_size = 0;
                    ifs.read((char *) &temp_size, sizeof(temp_size));
                    for (int source_node = node_begin; source_node < node_end; ++source_node) {
                        Grid_label &l = grid_labels[source_node];
                        int visible_size;
                        ifs.read((char *) &visible_size, sizeof(visible_size));
                        l.hub_labels.resize(visible_size);
                        ifs.read((char *) &l.a_p, sizeof(l.a_p));
                        ifs.read((char *) &l.b_p, sizeof(l.b_p));
                        ifs.read((char *) &l.c_p, sizeof(l.c_p));
                        ifs.read((char *) &l.d_p, sizeof(l.d_p));
                        if (visible_size > 0) {
                            for (int i = 0; i < visible_size; ++i) {
                                ifs.read((char *) &l.hub_labels[i].hub_id, sizeof(l.hub_labels[i].hub_id));
                                int convex_size;
                                ifs.read((char *) &convex_size, sizeof(convex_size));
                                l.hub_labels[i].convex_labels.resize(convex_size);
                                for (int j = 0; j < convex_size; ++j) {
                                    ifs.read((char *) &l.hub_labels[i].convex_labels[j].convex_vertex,
                                             sizeof(l.hub_labels[i].convex_labels[j].convex_vertex));
                                    ifs.read((char *) &l.hub_labels[i].convex_labels[j].distance,
                                             sizeof(l.hub_labels[i].convex_labels[j].distance));
                                    ifs.read((char *) &l.hub_labels[i].convex_labels[j].visibility,
                                             sizeof(l.hub_labels[i].convex_labels[j].visibility));
                                    ifs.read((char *) &l.hub_labels[i].convex_labels[j].predecessor,
                                             sizeof(l.hub_labels[i].convex_labels[j].predecessor));
                                }
                            }
                        }
                    }
                    ifs.close();
                }
            }
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


        const vector<Point> &get_grid_corner_point(int grid_index) {
//            int y = floor(grid_index/num_of_columns);
//            int x = grid_index - y * num_of_columns;
//
            int x = floor(grid_index / num_of_rows);
            int y = grid_index - num_of_rows * x;
            //CCW
            grid_corner_points.push_back(Point{(double) x, (double) y});
            grid_corner_points.push_back(Point{(double) x + grid_size, (double) y});
            grid_corner_points.push_back(Point{(double) x + grid_size, (double) y + grid_size});
            grid_corner_points.push_back(Point{(double) x, (double) y + grid_size});
            return grid_corner_points;
        }

        int  get_grid_index(Point p) {
            return floor(p.x/grid_size) * num_of_rows + floor(p.y/grid_size);
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

    };


}




#endif //LEARNED_HEURISTIC_EBHL_H
