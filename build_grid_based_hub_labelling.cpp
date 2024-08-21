/*
 * EHL construction and implementation
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <expansion.h>
#include "visibleAreaSearchInstance.h"
#include "point.h"
#include <omp.h>
#include <iomanip>
#include "geometry.h"
#include "ebhl.h"
#include <map>
#include <set>
#include <searchinstance.h>
#include <helper.h>
#include "coverage_ordering_path.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>


typedef boost::geometry::model::d2::point_xy<double> point_type;
typedef boost::geometry::model::polygon<point_type> polygon_type;
typedef boost::geometry::model::box<point_type> box_type;
using namespace std;
namespace pl = polyanya;
pl::MeshPtr mp;
pl::visibleAreaSearchInstance* vs;
pl::SearchInstance* se;
namespace bg = boost::geometry;

struct vis_poly {
    polygon_type  boost_poly;
    box_type bounding_box;
};


//output triangles into an entire polygon for visible area
void output_non_taut_triangle(const vector<pair<vis_poly,vis_poly>> &boost_poly_set, int &size, string output_file){
    cout<< "Saving triangles file to " << output_file << endl;
    ofstream outputFile(output_file);
    outputFile<<boost_poly_set.size()<<std::endl;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < boost_poly_set[i].first.boost_poly.outer().size(); ++j){
            outputFile<<std::fixed<<setprecision(8) <<i*2<<" "<<boost_poly_set[i].first.boost_poly.outer()[j].x()<<" "<<boost_poly_set[i].first.boost_poly.outer()[j].y() <<"\n";
        }
        for(int j = 0; j < boost_poly_set[i].second.boost_poly.outer().size(); ++j){
            outputFile<<std::fixed<<setprecision(8) <<i*2+1<<" "<<boost_poly_set[i].second.boost_poly.outer()[j].x()<<" "<<boost_poly_set[i].second.boost_poly.outer()[j].y() <<"\n";
        }
    }
    outputFile.close();

}

//checks whether a grid is within a rectangle or not
bool is_within_rectangle( const point_type& check_point, const point_type&  min, const point_type& max){

    return min.x()<= check_point.x() && check_point.x() <= max.x() && min.y()<= check_point.y() && check_point.y() <= max.y();

}

//checks whether a polygon is within a grid
bool polygon_within_grid(const polygon_type& poly1, const polygon_type& poly2){
    for(const auto& p : poly1.outer()){
        if(!is_within_rectangle(p,poly2.outer()[0],poly2.outer()[2] )){
            return false;
        }
    }
    return true;
}




bool are_rectangles_overlap( const point_type& r1_min,const point_type& r1_max, const point_type&  r2_min, const point_type& r2_max){

    bool overlap_case = r1_min.x() <= r2_max.x() && r2_min.x() <= r1_max.x() &&
                        r1_min.y()<= r2_max.y()&& r2_min.y() <= r1_max.y();

    bool contain_case1 = r1_min.x() >= r2_min.x() && r1_max.x() <= r2_max.x() &&
                         r1_min.y() >= r2_min.y() && r1_max.y() <= r2_max.y();

    bool contain_case2 = r2_min.x() >= r1_min.x() && r2_max.x() <= r1_max.x() &&
                         r2_min.y() >= r1_min.y() && r2_max.y() <= r1_max.y();

    return overlap_case || contain_case1 || contain_case2;

}

bool grid_within_polygon(const polygon_type& grid, const polygon_type& poly2){
    polygon_type max_grid;
    bg::append(max_grid,point_type{grid.outer()[0].x()-EPSILON, grid.outer()[0].y()-EPSILON});
    bg::append(max_grid,point_type{grid.outer()[1].x()-EPSILON, grid.outer()[1].y()+EPSILON});
    bg::append(max_grid,point_type{grid.outer()[2].x()+EPSILON, grid.outer()[2].y()+EPSILON});
    bg::append(max_grid,point_type{grid.outer()[3].x()+EPSILON, grid.outer()[3].y()-EPSILON});
    bg::append(max_grid,point_type{grid.outer()[4].x()-EPSILON, grid.outer()[4].y()-EPSILON});

    polygon_type min_grid;
    bg::append(min_grid,point_type{grid.outer()[0].x()+EPSILON, grid.outer()[0].y()+EPSILON});
    bg::append(min_grid,point_type{grid.outer()[1].x()+EPSILON, grid.outer()[1].y()-EPSILON});
    bg::append(min_grid,point_type{grid.outer()[2].x()-EPSILON, grid.outer()[2].y()-EPSILON});
    bg::append(min_grid,point_type{grid.outer()[3].x()-EPSILON, grid.outer()[3].y()+EPSILON});
    bg::append(min_grid,point_type{grid.outer()[4].x()+EPSILON, grid.outer()[4].y()+EPSILON});


    return   bg::within(max_grid,poly2) &&  bg::within(min_grid,poly2) &&  bg::within(grid,poly2);

    ;

}


//checks whether a given path is a non taut path could be pruned or not
bool is_non_taut_path( pl::Point start ,const pl::Vertex& vertex,const pl::Point& target, int is_left_bottom){
    // for testing grid only
    const vector<pl::Vertex>& mesh_vertices  = mp->mesh_vertices;
    const pl::Orientation& v_right_ori = get_orientation(start,vertex.p, target);
    if(v_right_ori == pl::Orientation::COLLINEAR){
        // collinear is an ambiguious case, set it to true;
        if(is_left_bottom){
            return !(start.distance(target) >= vertex.p.distance(target)
                     && start.distance(target) >= start.distance(vertex.p));
        }else{
            return true;
        }
    }
    pl::Point mid_point = {(mesh_vertices[vertex.obstacle_edge[0]].p.x+ mesh_vertices[vertex.obstacle_edge[1]].p.x)/2,
                           (mesh_vertices[vertex.obstacle_edge[0]].p.y+ mesh_vertices[vertex.obstacle_edge[1]].p.y)/2
    };
    const pl::Orientation& start_v_ori = get_orientation(start, vertex.p, mid_point);

    if(v_right_ori == start_v_ori){
        //same side;
        const pl::Orientation& o1 = get_orientation(mid_point, vertex.p, start);
        const pl::Orientation& o2 = get_orientation(mid_point, vertex.p, target);
        return o1 == o2;
    }else{
        return true;
    }
}

bool grid_to_label_is_non_taut(const pl::Point&  predecessor, const pl::Vertex&  convex_vertex, const pl::Grid_label& grid){

    return is_non_taut_path(predecessor, convex_vertex, grid.a_p,true) && is_non_taut_path(predecessor, convex_vertex, grid.b_p,false) &&
           is_non_taut_path(predecessor, convex_vertex, grid.c_p,false) && is_non_taut_path(predecessor, convex_vertex, grid.d_p,false);

}


bool is_turn_edge(int source_vertex, int target_vertex ){
    pl::Point start = mp->mesh_vertices[source_vertex].p;
    const pl::Vertex& v  = mp->mesh_vertices[target_vertex];
    const pl::Vertex& v1  = mp->mesh_vertices[v.obstacle_edge[0]];
    const pl::Vertex& v2  = mp->mesh_vertices[v.obstacle_edge[1]];
    const pl::Orientation  o1 = get_orientation(start,v.p,v1.p);
    const pl::Orientation  o2 = get_orientation(start,v.p,v2.p);
    if(o1==pl::Orientation::COLLINEAR || o2==pl::Orientation::COLLINEAR ){
        return true;
    }else if(o1 == o2 ){
        return true;
    }
    return false;
}

bool is_valid_turn_edge(int source_vertex, int target_vertex ){
    if(!mp->mesh_vertices[target_vertex].is_ambig){
        if(!is_turn_edge(source_vertex,target_vertex)){
            return false;
        }
    }

    if(!mp->mesh_vertices[source_vertex].is_ambig){
        if(!is_turn_edge(target_vertex,source_vertex)){
            return false;
        }
    }

    return true;

}


/*
 * Construct of EHL as detailed in Offline Preprocessing 3-5.
 * Using pre-constructed visibility graph and hub labels, EHL is built.
 *
 */
void build_ebhl(string dir, string map, int grid_size){

    std::cout<<"Building EHL ..."<< map << endl;
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

    //loading hub label
    PLabel lab;
    string p_label ="dataset/hub_label/" + dir + "/" + map + ".label";
    lab.load_labels(p_label.c_str());
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ifstream order_ifs(order.c_str());
    vector<NodeID> label_rank(numOfVertices);
    vector<NodeID> label_inv(numOfVertices);
    for(int i = 0; i < numOfVertices; ++i){
        NodeID tv;
        order_ifs >> tv;
        label_rank[tv] = i;
        label_inv[i] = tv;
    }

    //Pruning 4: filtering dead-end label;
    vector<vector<raw_label>> r_label = lab.get_raw_label_list(label_inv);
    unsigned label_size  = 0;
    for(auto fl : r_label){
        label_size += fl.size();
    }
    std::cout<<"Before filtering label size: " << label_size <<std::endl;

    vector<vector<raw_label>> f_label;
    for(unsigned i = 0; i < r_label.size(); i ++){
        int current_vertex_id = turning_vertices[i];
        vector<raw_label> filtered_labels;
        for(unsigned j = 0 ; j < r_label[i].size(); j ++){
            int current_hub_id = turning_vertices[label_inv[r_label[i][j].hub_id]];
            int current_predecessor  = turning_vertices[r_label[i][j].label_predecessor];
            int current_first_node  = turning_vertices[r_label[i][j].label_first_node];
            if(current_hub_id == current_vertex_id){
                filtered_labels.push_back(r_label[i][j]);
            }else{
                if(current_predecessor != current_vertex_id){
                    if(!is_valid_turn_edge(current_predecessor, current_vertex_id)){
                        continue;
                    }
                }else{
                    std::cout<< "current predecessor should not equal to current vertex id "<< std::endl;
                }
                if(current_first_node != current_hub_id ){
                    if(!is_valid_turn_edge(current_first_node,current_hub_id)){
                        continue;
                    }
                }else{
                    std::cout<< "current first node should not equal to current vertex id "<< std::endl;
                }
                filtered_labels.push_back(r_label[i][j]);
            }
        }
        filtered_labels.shrink_to_fit();
        f_label.push_back(filtered_labels);
    }
    f_label.shrink_to_fit();
    label_size  = 0;
    for(auto fl : f_label){
        label_size += fl.size();
    }
    std::cout<<"After filtering label size: " << label_size <<std::endl;

    int map_height;
    int map_width;
    mp->get_grid_width_height(map_width, map_height);
    //superimpose uniform grid for the whole map given the map size
    ebhl->initialize_grid_map(grid_size, map_width, map_height);

    //find all the triangle successors and insert into a unique set for computation of polygon later
    vector<pl::visibleAreaSearchInstance* > thread_search(omp_get_max_threads());
    for(unsigned i = 0; i < thread_search.size(); i ++ ){
        thread_search[i] = new pl::visibleAreaSearchInstance(mp);
    }
    vector< pair <vis_poly , vis_poly>> boost_poly (turning_point.size());
    //finding all the visible areas for the given map and pruning 1 is applied
    int progress = 0;
    {
        printf("Using %d threads\n", omp_get_max_threads());
#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
            const int node_count = turning_point.size();
            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            auto search = thread_search[thread_id];
            for(int source_node=node_begin; source_node < node_end; ++source_node) {
                search->search_non_taut_visible_area(turning_vertices[source_node], turning_vertices_location[source_node],boost_poly[source_node].first.boost_poly,boost_poly[source_node].second.boost_poly);
                //finds the two visible polygon for each convex vertex and stores it as two boost polygons
                bg::envelope(boost_poly[source_node].first.boost_poly,  boost_poly[source_node].first.bounding_box);
                bg::dsv( boost_poly[source_node].first.bounding_box);
                bg::envelope(boost_poly[source_node].second.boost_poly,  boost_poly[source_node].second.bounding_box);
                bg::dsv( boost_poly[source_node].second.bounding_box);
#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / node_count * 100.0;
                        std::cout << "Progress: [" << progress << "/" << node_count  << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }
            }
        }
    }
    //find the hub nodes and via labels for each grid
    progress = 0;
    {
        printf("Using %d threads\n", omp_get_max_threads());
#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
            const int node_count = ebhl->grid_labels.size();
            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            for(int source_node=node_begin; source_node < node_end; ++source_node) {
                // compute visible and partial visible vertices for each grid.
                vector<int> visible_vertices;
                vector<int> partial_visible;
                polygon_type grid;
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].a_p.x, ebhl->grid_labels[source_node].a_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].d_p.x, ebhl->grid_labels[source_node].d_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].c_p.x, ebhl->grid_labels[source_node].c_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].b_p.x, ebhl->grid_labels[source_node].b_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].a_p.x, ebhl->grid_labels[source_node].a_p.y));

                const point_type& grid_min = grid.outer()[0];
                const point_type& grid_max = grid.outer()[2];
                //checks whether it is partially or fully visible for first visible polygon
                for (unsigned j = 0; j <  boost_poly.size(); j++) {
                    bool already_added_vertex= false;
                    const point_type& min = boost_poly[j].first.bounding_box.min_corner();
                    const point_type& max = boost_poly[j].first.bounding_box.max_corner();

                    if(are_rectangles_overlap(min,max,grid_min,grid_max) ) {
                        //Fixed parameters passed incorrectly
                        if (polygon_within_grid(boost_poly[j].first.boost_poly,grid)) {
                            partial_visible.push_back(j);
                            already_added_vertex = true;
                        } else {
                            if (bg::intersects(grid, boost_poly[j].first.boost_poly)) {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].first.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (bg::covered_by(grid, boost_poly[j].first.boost_poly) && all_points_inside_grids) {
                                    visible_vertices.push_back(j);
                                }else{
                                    partial_visible.push_back(j);
                                }
                                already_added_vertex = true;
                            } else {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].first.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (all_points_inside_grids) {
                                    visible_vertices.push_back(j);
                                    already_added_vertex = true;

                                }

                            }

                        }
                    }
                    if(already_added_vertex){ continue; }

                    const point_type& second_min = boost_poly[j].second.bounding_box.min_corner();
                    const point_type& second_max = boost_poly[j].second.bounding_box.max_corner();
                    //checks whether its visible for second polygon
                    if(are_rectangles_overlap(second_min ,second_max,grid_min,grid_max)) {
                        //Fixed parameters passed incorrectly
                        if (polygon_within_grid(boost_poly[j].second.boost_poly,grid)) {
                            partial_visible.push_back(j);
                        } else {
                            if (bg::intersects(grid, boost_poly[j].second.boost_poly)) {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].second.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (bg::covered_by(grid, boost_poly[j].second.boost_poly) && all_points_inside_grids ) {
                                    visible_vertices.push_back(j);

                                }else{
                                    partial_visible.push_back(j);
                                }

                            } else {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].second.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (all_points_inside_grids) {
                                    visible_vertices.push_back(j);
                                }
                            }
                        }
                    }
                }

                visible_vertices.shrink_to_fit();
                partial_visible.shrink_to_fit();


                std::map<int, polyanya::Hub_label> label_mapper;
                //organises and inserts visible via labels for each hub node
                if(!visible_vertices.empty()){
                    for(unsigned j = 0 ; j < visible_vertices.size(); j ++){
                        const vector<raw_label>& label_list  = f_label[visible_vertices[j]];
                        for(unsigned k = 0 ; k < label_list.size(); k ++){
                            const int& hub_id = label_list[k].hub_id;
                            polyanya::Convex_vertices_label convex_label = polyanya::Convex_vertices_label{visible_vertices[j],label_list[k].label_distance,true, label_list[k].label_predecessor};
                            if(label_mapper.find(hub_id) != label_mapper.end()){
                                label_mapper[hub_id ].convex_labels.push_back(convex_label);
                            }else{
                                label_mapper.insert({hub_id, pl::Hub_label{hub_id,vector<pl::Convex_vertices_label>(0)} });
                                label_mapper[hub_id].convex_labels.push_back(convex_label);
                            }
                        }
                    }
                }

                //organises and inserts partially visible via labels for each hub node
                if(!partial_visible.empty()){
                    for(unsigned j = 0 ; j < partial_visible.size(); j ++){
                        const vector<raw_label>& label_list  = f_label[partial_visible[j]];
                        for(unsigned k = 0 ; k < label_list.size(); k ++){
                            const int& hub_id = label_list[k].hub_id;
                            polyanya::Convex_vertices_label convex_label = polyanya::Convex_vertices_label{partial_visible[j],label_list[k].label_distance,false, label_list[k].label_predecessor};
                            if(label_mapper.find(hub_id) != label_mapper.end()){
                                label_mapper[hub_id ].convex_labels.push_back(convex_label);
                            }else{
                                label_mapper.insert({hub_id, pl::Hub_label{hub_id,vector<pl::Convex_vertices_label>(0)} });
                                label_mapper[hub_id].convex_labels.push_back(convex_label);
                            }
                        }
                    }
                }

                ebhl->grid_labels[source_node].hub_labels = vector<pl::Hub_label>() ;

                std::map<int,pl::Hub_label>::iterator it;

                for (it = label_mapper.begin(); it != label_mapper.end(); it++)
                {
                    ebhl->grid_labels[source_node].hub_labels.push_back(it->second);
                }
                //sort the labels based on hub node ID
                std::sort(std::begin(ebhl->grid_labels[source_node].hub_labels),
                          std::end(ebhl->grid_labels[source_node].hub_labels),
                          []( const pl::Hub_label&  lvalue, const pl::Hub_label& rvalue) {
                              return  lvalue.hub_id < rvalue.hub_id; });

                const pl::Point &  min_corner =  ebhl->grid_labels[source_node].a_p;
                const pl::Point &  max_corner =  ebhl->grid_labels[source_node].c_p;
                //Pruning 2 and 3
                for( unsigned j = 0;  j < ebhl->grid_labels[source_node].hub_labels.size(); j ++ ){
                    const pl::Hub_label hl = ebhl->grid_labels[source_node].hub_labels[j];
                    const vector< pl::Convex_vertices_label>& convex_label_list = hl.convex_labels;
                    vector< pl::Convex_vertices_label> filtered_convex_label_list;
                    double min_upperbound = INF_WEIGHT ;
                    for(const pl::Convex_vertices_label & label : convex_label_list){
                        if(label.visibility ){
                            // distance function should already handle the containment cases.
                            double upperbound = turning_point[label.convex_vertex].max_distance_to_grid( min_corner, max_corner) + label.distance;
                            if(min_upperbound > upperbound){
                                min_upperbound = upperbound;
                            }
                        }
                    }

                    for(const pl::Convex_vertices_label& label : convex_label_list){
                        double lowerbound = turning_point[label.convex_vertex].min_distance_to_grid( min_corner, max_corner) + label.distance;
                        if(lowerbound <= min_upperbound){
                            //only add label is  lowerbound  <= min_upperbound
                            if(label.visibility){
                                // prune non-taut path label ;
                                pl::Point hub =  turning_point[label_inv[hl.hub_id]];
                                pl::Point predecessor =  turning_point[label.predecessor];
                                const pl::Vertex& convex_vertex = mp->mesh_vertices[turning_vertices[label.convex_vertex]];
                                if(hub == convex_vertex.p){
                                    filtered_convex_label_list.push_back(label);
                                }else{
                                    if(!grid_to_label_is_non_taut(predecessor,convex_vertex, ebhl->grid_labels[source_node])) {
                                        filtered_convex_label_list.push_back(label);
                                    }
                                }
                            }else{
                                filtered_convex_label_list.push_back(label);
                            }
                        }
                    }
                    ebhl->grid_labels[source_node].hub_labels[j].convex_labels = filtered_convex_label_list;
                    ebhl->grid_labels[source_node].hub_labels[j].convex_labels.shrink_to_fit();
                }
                ebhl->grid_labels[source_node].hub_labels.erase(std::remove_if( ebhl->grid_labels[source_node].hub_labels.begin(), ebhl->grid_labels[source_node].hub_labels.end(),
                                                                      [](const pl::Hub_label& x) {
                                                                          return x.convex_labels.empty(); // put your condition here
                                                                      }), ebhl->grid_labels[source_node].hub_labels.end());
                pl::Convex_vertices_label cv = pl::Convex_vertices_label{ -1,INF_WEIGHT,false, -1};
                polyanya::Hub_label hub_label = polyanya::Hub_label{ (int)turning_point.size(), vector<pl::Convex_vertices_label>{cv}};
                ebhl->grid_labels[source_node].hub_labels.push_back( hub_label );

                //Optimisation: storing the lower bound;
                for(unsigned j = 0; j <ebhl->grid_labels[source_node].hub_labels.size(); j ++){
                    double lower_bound = std::numeric_limits<double>::max();
                    if(ebhl->grid_labels[source_node].hub_labels[j].hub_id ==(int)turning_point.size() ){
                        ebhl->grid_labels[source_node].hub_labels[j].min_lower_bound = lower_bound;
                        continue;
                    }
                    for(auto convex_label : ebhl->grid_labels[source_node].hub_labels[j].convex_labels ){
                        lower_bound =  std::min(lower_bound,turning_point[convex_label.convex_vertex].min_distance_to_grid( min_corner, max_corner) + convex_label.distance);
                    }
                    ebhl->grid_labels[source_node].hub_labels[j].min_lower_bound = lower_bound;
                }
                ebhl->grid_labels[source_node].hub_labels.shrink_to_fit();

#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / node_count * 100.0;
                        std::cout << "Progress: [" << progress << "/" << node_count  << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }

            }
        }
    }

    timer.stop();
    std::cout<<std::fixed << setprecision(8) <<"Preprocessing finished in "<< timer.elapsed_time_micro()/1000000 << " seconds"<< std::endl;

    //save in adjacent lists without parallelisation
    string output_path = "dataset/ehl/" + dir + "/" + map + ".nt_adj_"+ to_string(grid_size);
    ebhl->save_adjacent_list(output_path.c_str());
    //store the visible polygons
    string output_triangle_file = "dataset/ehl/" + dir + "/" + map +".nt_triangles_"+ to_string(grid_size);
    int size_of_turning = turning_point.size();
    output_non_taut_triangle(boost_poly, size_of_turning,output_triangle_file);

    delete ebhl;

}

/*
 * Construct of EHL as detailed in Offline Preprocessing 3-5.
 * Using pre-constructed visibility graph and hub labels, EHL is built.
 * Same construction as EHL, except saving it as grid labels (not adjacent list) for ease of loading in building EHL*
 *
 */
void build_ebhl_workload(string dir, string map, int grid_size){

    std::cout<<"Building EHL ..."<< map << endl;
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

    //loading hub label
    PLabel lab;
    string p_label ="dataset/hub_label/" + dir + "/" + map + ".label";
    lab.load_labels(p_label.c_str());
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ifstream order_ifs(order.c_str());
    vector<NodeID> label_rank(numOfVertices);
    vector<NodeID> label_inv(numOfVertices);
    for(int i = 0; i < numOfVertices; ++i){
        NodeID tv;
        order_ifs >> tv;
        label_rank[tv] = i;
        label_inv[i] = tv;
    }

    //Pruning 4: filtering dead-end label;
    vector<vector<raw_label>> r_label = lab.get_raw_label_list(label_inv);
    unsigned label_size  = 0;
    for(auto fl : r_label){
        label_size += fl.size();
    }
    std::cout<<"Before filtering label size: " << label_size <<std::endl;

    vector<vector<raw_label>> f_label;
    for(unsigned i = 0; i < r_label.size(); i ++){
        int current_vertex_id = turning_vertices[i];
        vector<raw_label> filtered_labels;
        for(unsigned j = 0 ; j < r_label[i].size(); j ++){
            int current_hub_id = turning_vertices[label_inv[r_label[i][j].hub_id]];
            int current_predecessor  = turning_vertices[r_label[i][j].label_predecessor];
            int current_first_node  = turning_vertices[r_label[i][j].label_first_node];
            if(current_hub_id == current_vertex_id){
                filtered_labels.push_back(r_label[i][j]);
            }else{
                if(current_predecessor != current_vertex_id){
                    if(!is_valid_turn_edge(current_predecessor, current_vertex_id)){
                        continue;
                    }
                }else{
                    std::cout<< "current predecessor should not equal to current vertex id "<< std::endl;
                }
                if(current_first_node != current_hub_id ){
                    if(!is_valid_turn_edge(current_first_node,current_hub_id)){
                        continue;
                    }
                }else{
                    std::cout<< "current first node should not equal to current vertex id "<< std::endl;
                }
                filtered_labels.push_back(r_label[i][j]);
            }
        }
        filtered_labels.shrink_to_fit();
        f_label.push_back(filtered_labels);
    }
    f_label.shrink_to_fit();
    label_size  = 0;
    for(auto fl : f_label){
        label_size += fl.size();
    }
    std::cout<<"After filtering label size: " << label_size <<std::endl;

    int map_height;
    int map_width;
    mp->get_grid_width_height(map_width, map_height);
    //superimpose uniform grid for the whole map given the map size
    ebhl->initialize_grid_map(grid_size, map_width, map_height);

    //find all the triangle successors and insert into a unique set for computation of polygon later
    vector<pl::visibleAreaSearchInstance* > thread_search(omp_get_max_threads());
    for(unsigned i = 0; i < thread_search.size(); i ++ ){
        thread_search[i] = new pl::visibleAreaSearchInstance(mp);
    }
    vector< pair <vis_poly , vis_poly>> boost_poly (turning_point.size());
    //finding all the visible areas for the given map and pruning 1 is applied
    int progress = 0;
    {
        printf("Using %d threads\n", omp_get_max_threads());
#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
            const int node_count = turning_point.size();
            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            auto search = thread_search[thread_id];
            for(int source_node=node_begin; source_node < node_end; ++source_node) {
                search->search_non_taut_visible_area(turning_vertices[source_node], turning_vertices_location[source_node],boost_poly[source_node].first.boost_poly,boost_poly[source_node].second.boost_poly);
                //finds the two visible polygon for each convex vertex and stores it as two boost polygons
                bg::envelope(boost_poly[source_node].first.boost_poly,  boost_poly[source_node].first.bounding_box);
                bg::dsv( boost_poly[source_node].first.bounding_box);
                bg::envelope(boost_poly[source_node].second.boost_poly,  boost_poly[source_node].second.bounding_box);
                bg::dsv( boost_poly[source_node].second.bounding_box);
#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / node_count * 100.0;
                        std::cout << "Progress: [" << progress << "/" << node_count  << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }
            }
        }
    }
    //find the hub nodes and via labels for each grid
    progress = 0;
    {
        printf("Using %d threads\n", omp_get_max_threads());
#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
            const int node_count = ebhl->grid_labels.size();
            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            for(int source_node=node_begin; source_node < node_end; ++source_node) {
                // compute visible and partial visible vertices for each grid.
                vector<int> visible_vertices;
                vector<int> partial_visible;
                polygon_type grid;
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].a_p.x, ebhl->grid_labels[source_node].a_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].d_p.x, ebhl->grid_labels[source_node].d_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].c_p.x, ebhl->grid_labels[source_node].c_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].b_p.x, ebhl->grid_labels[source_node].b_p.y));
                bg::append(grid.outer(), point_type(ebhl->grid_labels[source_node].a_p.x, ebhl->grid_labels[source_node].a_p.y));

                const point_type& grid_min = grid.outer()[0];
                const point_type& grid_max = grid.outer()[2];
                //checks whether it is partially or fully visible for first visible polygon
                for (unsigned j = 0; j <  boost_poly.size(); j++) {
                    bool already_added_vertex= false;
                    const point_type& min = boost_poly[j].first.bounding_box.min_corner();
                    const point_type& max = boost_poly[j].first.bounding_box.max_corner();

                    if(are_rectangles_overlap(min,max,grid_min,grid_max) ) {
                        //Fixed parameters passed incorrectly
                        if (polygon_within_grid(boost_poly[j].first.boost_poly,grid)) {
                            partial_visible.push_back(j);
                            already_added_vertex = true;
                        } else {
                            if (bg::intersects(grid, boost_poly[j].first.boost_poly)) {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].first.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (bg::covered_by(grid, boost_poly[j].first.boost_poly) && all_points_inside_grids) {
                                    visible_vertices.push_back(j);
                                }else{
                                    partial_visible.push_back(j);
                                }
                                already_added_vertex = true;
                            } else {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].first.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (all_points_inside_grids) {
                                    visible_vertices.push_back(j);
                                    already_added_vertex = true;

                                }

                            }

                        }
                    }
                    if(already_added_vertex){ continue; }

                    const point_type& second_min = boost_poly[j].second.bounding_box.min_corner();
                    const point_type& second_max = boost_poly[j].second.bounding_box.max_corner();
                    //checks whether its visible for second polygon
                    if(are_rectangles_overlap(second_min ,second_max,grid_min,grid_max)) {
                        //Fixed parameters passed incorrectly
                        if (polygon_within_grid(boost_poly[j].second.boost_poly,grid)) {
                            partial_visible.push_back(j);
                        } else {
                            if (bg::intersects(grid, boost_poly[j].second.boost_poly)) {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].second.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (bg::covered_by(grid, boost_poly[j].second.boost_poly) && all_points_inside_grids ) {
                                    visible_vertices.push_back(j);

                                }else{
                                    partial_visible.push_back(j);
                                }

                            } else {
                                bool all_points_inside_grids = true;
                                for(const auto& gp : grid.outer()){
                                    if(!bg::covered_by(gp,boost_poly[j].second.boost_poly)){
                                        all_points_inside_grids = false;
                                        break;
                                    }
                                }
                                if (all_points_inside_grids) {
                                    visible_vertices.push_back(j);
                                }
                            }
                        }
                    }
                }

                visible_vertices.shrink_to_fit();
                partial_visible.shrink_to_fit();


                std::map<int, polyanya::Hub_label> label_mapper;
                //organises and inserts visible via labels for each hub node
                if(!visible_vertices.empty()){
                    for(unsigned j = 0 ; j < visible_vertices.size(); j ++){
                        const vector<raw_label>& label_list  = f_label[visible_vertices[j]];
                        for(unsigned k = 0 ; k < label_list.size(); k ++){
                            const int& hub_id = label_list[k].hub_id;
                            polyanya::Convex_vertices_label convex_label = polyanya::Convex_vertices_label{visible_vertices[j],label_list[k].label_distance,true, label_list[k].label_predecessor};
                            if(label_mapper.find(hub_id) != label_mapper.end()){
                                label_mapper[hub_id ].convex_labels.push_back(convex_label);
                            }else{
                                label_mapper.insert({hub_id, pl::Hub_label{hub_id,vector<pl::Convex_vertices_label>(0)} });
                                label_mapper[hub_id].convex_labels.push_back(convex_label);
                            }
                        }
                    }
                }

                //organises and inserts partially visible via labels for each hub node
                if(!partial_visible.empty()){
                    for(unsigned j = 0 ; j < partial_visible.size(); j ++){
                        const vector<raw_label>& label_list  = f_label[partial_visible[j]];
                        for(unsigned k = 0 ; k < label_list.size(); k ++){
                            const int& hub_id = label_list[k].hub_id;
                            polyanya::Convex_vertices_label convex_label = polyanya::Convex_vertices_label{partial_visible[j],label_list[k].label_distance,false, label_list[k].label_predecessor};
                            if(label_mapper.find(hub_id) != label_mapper.end()){
                                label_mapper[hub_id ].convex_labels.push_back(convex_label);
                            }else{
                                label_mapper.insert({hub_id, pl::Hub_label{hub_id,vector<pl::Convex_vertices_label>(0)} });
                                label_mapper[hub_id].convex_labels.push_back(convex_label);
                            }
                        }
                    }
                }

                ebhl->grid_labels[source_node].hub_labels = vector<pl::Hub_label>() ;

                std::map<int,pl::Hub_label>::iterator it;

                for (it = label_mapper.begin(); it != label_mapper.end(); it++)
                {
                    ebhl->grid_labels[source_node].hub_labels.push_back(it->second);
                }
                //sort the labels based on hub node ID
                std::sort(std::begin(ebhl->grid_labels[source_node].hub_labels),
                          std::end(ebhl->grid_labels[source_node].hub_labels),
                          []( const pl::Hub_label&  lvalue, const pl::Hub_label& rvalue) {
                              return  lvalue.hub_id < rvalue.hub_id; });

                const pl::Point &  min_corner =  ebhl->grid_labels[source_node].a_p;
                const pl::Point &  max_corner =  ebhl->grid_labels[source_node].c_p;
                //Pruning 2 and 3
                for( unsigned j = 0;  j < ebhl->grid_labels[source_node].hub_labels.size(); j ++ ){
                    const pl::Hub_label hl = ebhl->grid_labels[source_node].hub_labels[j];
                    const vector< pl::Convex_vertices_label>& convex_label_list = hl.convex_labels;
                    vector< pl::Convex_vertices_label> filtered_convex_label_list;
                    double min_upperbound = INF_WEIGHT ;
                    for(const pl::Convex_vertices_label & label : convex_label_list){
                        if(label.visibility ){
                            // distance function should already handle the containment cases.
                            double upperbound = turning_point[label.convex_vertex].max_distance_to_grid( min_corner, max_corner) + label.distance;
                            if(min_upperbound > upperbound){
                                min_upperbound = upperbound;
                            }
                        }
                    }

                    for(const pl::Convex_vertices_label& label : convex_label_list){
                        double lowerbound = turning_point[label.convex_vertex].min_distance_to_grid( min_corner, max_corner) + label.distance;
                        if(lowerbound <= min_upperbound){
                            //only add label is  lowerbound  <= min_upperbound
                            if(label.visibility){
                                // prune non-taut path label ;
                                pl::Point hub =  turning_point[label_inv[hl.hub_id]];
                                pl::Point predecessor =  turning_point[label.predecessor];
                                const pl::Vertex& convex_vertex = mp->mesh_vertices[turning_vertices[label.convex_vertex]];
                                if(hub == convex_vertex.p){
                                    filtered_convex_label_list.push_back(label);
                                }else{
                                    if(!grid_to_label_is_non_taut(predecessor,convex_vertex, ebhl->grid_labels[source_node])) {
                                        filtered_convex_label_list.push_back(label);
                                    }
                                }
                            }else{
                                filtered_convex_label_list.push_back(label);
                            }
                        }
                    }
                    ebhl->grid_labels[source_node].hub_labels[j].convex_labels = filtered_convex_label_list;
                    ebhl->grid_labels[source_node].hub_labels[j].convex_labels.shrink_to_fit();
                }
                ebhl->grid_labels[source_node].hub_labels.erase(std::remove_if( ebhl->grid_labels[source_node].hub_labels.begin(), ebhl->grid_labels[source_node].hub_labels.end(),
                                                                                [](const pl::Hub_label& x) {
                                                                                    return x.convex_labels.empty(); // put your condition here
                                                                                }), ebhl->grid_labels[source_node].hub_labels.end());
                pl::Convex_vertices_label cv = pl::Convex_vertices_label{ -1,INF_WEIGHT,false, -1};
                polyanya::Hub_label hub_label = polyanya::Hub_label{ (int)turning_point.size(), vector<pl::Convex_vertices_label>{cv}};
                ebhl->grid_labels[source_node].hub_labels.push_back( hub_label );

                //Optimisation: storing the lower bound;
                for(unsigned j = 0; j <ebhl->grid_labels[source_node].hub_labels.size(); j ++){
                    double lower_bound = std::numeric_limits<double>::max();
                    if(ebhl->grid_labels[source_node].hub_labels[j].hub_id ==(int)turning_point.size() ){
                        ebhl->grid_labels[source_node].hub_labels[j].min_lower_bound = lower_bound;
                        continue;
                    }
                    for(auto convex_label : ebhl->grid_labels[source_node].hub_labels[j].convex_labels ){
                        lower_bound =  std::min(lower_bound,turning_point[convex_label.convex_vertex].min_distance_to_grid( min_corner, max_corner) + convex_label.distance);
                    }
                    ebhl->grid_labels[source_node].hub_labels[j].min_lower_bound = lower_bound;
                }
                ebhl->grid_labels[source_node].hub_labels.shrink_to_fit();

#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / node_count * 100.0;
                        std::cout << "Progress: [" << progress << "/" << node_count  << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }

            }
        }
    }

    timer.stop();
    std::cout<<std::fixed << setprecision(8) <<"Preprocessing finished in "<< timer.elapsed_time_micro()/1000000 << " seconds"<< std::endl;

    //save in grid label without parallelisation
    string output_path = "dataset/ehl/" + dir + "/" + map + ".nt_gridlabel_"+ to_string(grid_size);
    ebhl->save_grid_labels(output_path.c_str());
    //store the visible polygons
    string output_triangle_file = "dataset/ehl/" + dir + "/" + map +".nt_triangles_"+ to_string(grid_size);
    int size_of_turning = turning_point.size();
    output_non_taut_triangle(boost_poly, size_of_turning,output_triangle_file);

    delete ebhl;

}


int main(int argc, char*argv[]){

    try{
        string dir;
        string map;
        string grid_size;
        if(argc != 4){
            cerr << argv[0] << "directory map grid_size" << endl;
            return 1;
        }else{
            dir = argv[1];
            map = argv[2];
            grid_size = argv[3];
        }

        //improved memory and preprocessing pruning version saving in adjacent lists
        //build_ebhl(dir,map,stoi(grid_size));

        //build EHL for EHL*
        build_ebhl_workload(dir,map,stoi(grid_size));


    }catch(exception&err){

        cerr << "Stopped on exception lol: " << err.what() << endl;
    }
}