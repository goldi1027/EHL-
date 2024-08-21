#include <stdio.h>
#include <iostream>
#include <fstream>
#include "visibleAreaSearchInstance.h"
#include <query_workload.h>
typedef boost::geometry::model::d2::point_xy<double> point_type;
typedef boost::geometry::model::polygon<point_type> polygon_type;
typedef boost::geometry::model::box<point_type> box_type;
using namespace std;
namespace pl = polyanya;
using namespace polyanya;
pl::MeshPtr mp;


/*
 * Output the mapper for EHL* into a text file
 */

void output_mapper(vector<int> mapper,string output_file) {

    // Open the binary file for writing
    std::ofstream outputFile(output_file);


    // Write the size of the vector to the file
    size_t size = mapper.size();
    outputFile.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // Write the vector data to the file
    outputFile.write(reinterpret_cast<const char*>(mapper.data()), size * sizeof(int));

    // Close the file
    outputFile.close();

}
/*
 * Constructs EHL* for clustered scenarios
 * Using pre-computed EHL, it loads the EHL index, calculates the memory and iteratively compresses it to the target memory
 */
void build_query_workload_index_load(string dir, string map, int grid_size, unsigned memory, string num_rects){
    std::cout<<"Building EHL Clustered..."<< map << endl;
    //load mesh
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    // mark obstacle edge;
    mp->pre_compute_obstacle_edge_on_vertex();
    // mark valid turning point;
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());
    //mp->mark_turning_point_polygon();
    //turning point is the actual point which would be used in search
    vector<pl::Point> turning_point;
    //corresponding vertices number in mesh_vertices
    vector<int> turning_vertices;
    //vertice location in polygon?
    vector<pl::PointLocation> turning_vertices_location;

    warthog::timer timer =  warthog::timer ();
    timer.start();
    int id = 0;
    int poly = 0;
    for( pl::Vertex& v : mp->mesh_vertices){
        if(v.is_turning_vertex && !v.is_ambig) {
            for (int polygon: v.polygons) {
                if (polygon != -1) {
//                    p.polygons.insert(polygon);
                    poly = polygon;
                }
            }

            // there is some issue here, old implementation assume that these vertex always inside one polygon only. it doesnt
            // use the correct type actually, I fix it here manually assign a random adjacent polygon
            pl::PointLocation location = {pl::PointLocation::ON_CORNER_VERTEX_UNAMBIG, poly, -1, id, -1};
            turning_vertices_location.push_back(location);
            turning_vertices.push_back(id);
            turning_point.push_back(mp->mesh_vertices[id].p);
        }
        id ++;
    }

    auto ebhl = new pl::EBHL(grid_size,mp->get_map_height(),mp->get_map_width());
    string grid_label_path = "dataset/ehl/" + dir + "/" + map + ".nt_gridlabel_"+ to_string(grid_size);
    ebhl->load_grid_labels(grid_label_path.c_str());

    cout << "Finish loading EHL" << endl;

    string triangles_path = "dataset/ehl/" + dir + "/" + map + ".nt_triangles_"+ to_string(grid_size);
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    int map_height;
    int map_width;
    mp->get_grid_width_height(map_width, map_height);

    auto workload = new pl::Workload(map_height, map_width,grid_size,memory,ebhl);

    //load clustered query file
    string query_file = "dataset/query_workload/" + dir + "/" + map + "_" + num_rects + "_query_workload.csv";

    //load the clustered historical queries on the map
    cout << "Load query" << endl;
    int highest_workload = 0;
    workload->load_query(query_file,highest_workload);

    //determines the traversable grids of the map to exclude merging it
    int traversable = 0;
    for(int i = 0; i < ebhl->grid_labels.size(); ++i){
        PointLocation a = mp->get_point_location(ebhl->grid_labels[i].a_p);
        PointLocation b = mp->get_point_location(ebhl->grid_labels[i].b_p);
        PointLocation c = mp->get_point_location(ebhl->grid_labels[i].c_p);
        PointLocation d = mp->get_point_location(ebhl->grid_labels[i].d_p);

        if(a.type == polyanya::PointLocation::NOT_ON_MESH and b.type == polyanya::PointLocation::NOT_ON_MESH and
           c.type == polyanya::PointLocation::NOT_ON_MESH and d.type== polyanya::PointLocation::NOT_ON_MESH){
            workload->traversable[i] = false;
        }
        else{
            traversable++;
        }
    }

    workload->init();

    delete ebhl;

    //main compression phase of EHL*
    workload->merge();

    timer.stop();
    string output_file = "dataset/ehl/" + dir + "/workload/" + map + "_" + to_string(memory) + "_" + num_rects+ ".workload_label";
    workload->save_adjacent_workload(output_file.c_str());

    cout << "================Done: " << output_file << "===================" <<endl;
    string mapper_out = "dataset/ehl/" + dir + "/workload/" + map + "_" + to_string(memory) + "_" + num_rects +  ".workload_mapper";
    output_mapper(workload->mapper,mapper_out);
    double mapper_size = (8 + (workload->mapper.size() * 4))/(double)1000000;

    unsigned long long hub_nodes = workload->get_hub_nodes();
    unsigned long long total_memory = workload->get_total_memory();


    double avg_grid = workload->get_avg_grid();


    cout << setprecision(8) << "Total memory: " << total_memory/(double)1000000 << endl;
    cout << setprecision(8) << "Processing time: " << timer.elapsed_time_micro()/1000000 << endl;
    string out_final = "dataset/result/preprocessing/" +  dir + "/" + dir + "_aaai_" + to_string(memory) + "_" + num_rects + ".csv";
    std::ofstream myFile_final(out_final, std::ofstream::out | std::ofstream::app);
    myFile_final << map << "," << timer.elapsed_time_micro()/1000000 << "," << total_memory/(double)1000000 << "," << mapper_size  << "," << avg_grid << "," << highest_workload<<
    "," << (hub_nodes)/(double)1000000 << "," << (total_memory/(double)1000000)+mapper_size+((hub_nodes)/(double)1000000) << "\n";
    myFile_final.close();

    delete workload;

}


void build_query_workload_uniform(string dir, string map, int grid_size, unsigned memory, string num_rects){
    std::cout<<"Building EHL-W Uniform ..."<< map << endl;
    //load mesh
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    // mark obstacle edge;
    mp->pre_compute_obstacle_edge_on_vertex();
    // mark valid turning point;
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());
    //mp->mark_turning_point_polygon();
    //turning point is the actual point which would be used in search
    vector<pl::Point> turning_point;
    //corresponding vertices number in mesh_vertices
    vector<int> turning_vertices;
    //vertice location in polygon?
    vector<pl::PointLocation> turning_vertices_location;

    warthog::timer timer =  warthog::timer ();
    timer.start();
    int id = 0;
    int poly = 0;
    for( pl::Vertex& v : mp->mesh_vertices){
        if(v.is_turning_vertex && !v.is_ambig) {
            for (int polygon: v.polygons) {
                if (polygon != -1) {
//                    p.polygons.insert(polygon);
                    poly = polygon;
                }
            }

            // there is some issue here, old implementation assume that these vertex always inside one polygon only. it doesnt
            // use the correct type actually, I fix it here manually assign a random adjacent polygon
            pl::PointLocation location = {pl::PointLocation::ON_CORNER_VERTEX_UNAMBIG, poly, -1, id, -1};
            turning_vertices_location.push_back(location);
            turning_vertices.push_back(id);
            turning_point.push_back(mp->mesh_vertices[id].p);
        }
        id ++;
    }

    auto ebhl = new pl::EBHL(grid_size,mp->get_map_height(),mp->get_map_width());
    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_gridlabel_"+ to_string(grid_size);
    ebhl->load_grid_labels(grid_label_path.c_str());

    cout << "Finish loading EHL" << endl;

    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_"+ to_string(grid_size);
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    int map_height;
    int map_width;
    mp->get_grid_width_height(map_width, map_height);


    auto workload = new pl::Workload(map_height, map_width,grid_size,memory,ebhl);

    int traversable = 0;
    for(int i = 0; i < ebhl->grid_labels.size(); ++i){
        PointLocation a = mp->get_point_location(ebhl->grid_labels[i].a_p);
        PointLocation b = mp->get_point_location(ebhl->grid_labels[i].b_p);
        PointLocation c = mp->get_point_location(ebhl->grid_labels[i].c_p);
        PointLocation d = mp->get_point_location(ebhl->grid_labels[i].d_p);

        if(a.type == polyanya::PointLocation::NOT_ON_MESH and b.type == polyanya::PointLocation::NOT_ON_MESH and
           c.type == polyanya::PointLocation::NOT_ON_MESH and d.type== polyanya::PointLocation::NOT_ON_MESH){
            workload->traversable[i] = false;
        }
        else{
            traversable++;
        }
    }
    workload->init();

    delete ebhl;


    cout << "Convex vertices: " << turning_point.size() << endl;
    workload->merge();

    timer.stop();
    string output_file = "dataset/ebhl/" + dir + "/workload/" + map + "_" + to_string(memory) + "_" + num_rects+ ".workload_label";
    //workload->save_adjacent_workload(output_file.c_str());

    cout << "================Done: " << output_file << "===================" <<endl;
    string mapper_out = "dataset/ebhl/" + dir + "/workload/" + map + "_" + to_string(memory) + "_" + num_rects +  ".workload_mapper";
    output_mapper(workload->mapper,mapper_out);
    double mapper_size = (8 + (workload->mapper.size() * 4))/(double)1000000;

    unsigned long long hub_nodes = workload->get_hub_nodes();
    unsigned long long total_memory = workload->get_total_memory();


    double avg_grid = workload->get_avg_grid();

    int highest_workload = 0;

//    cout << setprecision(8) << "Total memory: " << total_memory/(double)1000000 << endl;
//    cout << setprecision(8) << "Processing time: " << timer.elapsed_time_micro()/1000000 << endl;
//    string out = "dataset/result/debug/"  + map + "_stats_" + to_string(memory) + "_" + num_rects + ".csv";
//    std::ofstream myFile(out);
//    myFile << timer.elapsed_time_micro()/1000000 << "," << total_memory/(double)1000000 << "," << mapper_size  << "," << avg_grid<< "\n";
//    myFile << workload->grids.size() << "\n";
//    myFile << "ID,Workload,Merged Status,Mapped Where,Merged grids size,X,Y" << "\n";
//    for(int i = 0; i < workload->grids.size(); ++i){
//        if(workload->traversable[i]){
//            myFile << i << "," << workload->grids[i].workload << "," << workload->grid_merge_status[i] << "," << workload->mapper[workload->grids[i].id] << ","
//            << workload->grids[i].merged_grids.size() << "," << workload->grids[i].gridLabel.a_p.x << "," << workload->grids[i].gridLabel.a_p.y  <<"\n";
//        }
//    }
//    myFile.close();

    cout << setprecision(8) << "Total memory: " << total_memory/(double)1000000 << endl;
    cout << setprecision(8) << "Processing time: " << timer.elapsed_time_micro()/1000000 << endl;
    string out_final = "dataset/result/preprocessing/" +  dir + "/" + dir + "_aaai_" + to_string(memory) + "_" + num_rects + ".csv";
    std::ofstream myFile_final(out_final, std::ofstream::out | std::ofstream::app);
    myFile_final << map << "," << timer.elapsed_time_micro()/1000000 << "," << total_memory/(double)1000000 << "," << mapper_size  << "," << avg_grid << "," << highest_workload<<
                 "," << (hub_nodes)/(double)1000000 << "," << (total_memory/(double)1000000)+mapper_size+((hub_nodes)/(double)1000000) << "\n";
    myFile_final.close();

    delete workload;

}

int main(int argc, char*argv[]){

    try{
        string dir;
        string map;
        string grid_size;
        string memory;
        string num_rects;
        if(argc != 6){
            cerr << argv[0] << "directory map" << endl;
            return 1;
        }else{
            dir = argv[1];
            map = argv[2];
            grid_size = argv[3];
            memory = argv[4];
            num_rects = argv[5];
        }



        build_query_workload_index_load(dir,map,stoi(grid_size),stoi(memory),num_rects);

        //build_query_workload_uniform(dir,map,stoi(grid_size),stoi(memory),num_rects);

    }catch(exception&err){

        cerr << "Stopped on exception lol: " << err.what() << endl;
    }
}