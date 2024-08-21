/*
 * Implementation of query class for EHL
 * This class initialises and prepares each search given a start and target coordinate in the Euclidean space. For labels which are not visible, if needed, visibility test is also performed.
 */


#ifndef LEARNED_HEURISTIC_EBHL_QUERY_V2_H
#define LEARNED_HEURISTIC_EBHL_QUERY_V2_H

#include <mesh.h>
#include <ebhl.h>
#include "expansion.h"
#include "coverage_ordering_path.h"

namespace polyanya {
    typedef Mesh* MeshPtr;
    typedef EBHL* EBHLPtrV2;
// Polyanya instance for point to point search
    class EBHL_query_v2
    {
    private:
        const MeshPtr mesh;
        const EBHLPtrV2 ebhl;
        warthog::timer timer;

        std::vector<int> polygons_check_list;
        int visible_search_id = 0;
        std::vector<int> start_visibility_times_flag;
        std::vector<int> goal_visibility_times_flag;
        std::vector<bool> start_visibility_cache;
        std::vector<bool> goal_visibility_cache;
    public:

        Point start, goal;
        int start_index, goal_index;
        double distance;

        PLabel shp_lab;
        vector<int> shp_rank;
        vector<int> shp_inv;
        vector<Point> turning_point;
        vector<Point> obstacle_middle_point;

        int best_hub_id;
        int target_predecessor;
        int target_convex;

        int start_predecessor;
        int start_convex;

        vector<Point> shortest_path;

        EBHL_query_v2() = default;
        EBHL_query_v2(MeshPtr m , EBHLPtrV2 eb,  vector<Point> tp, vector<Point> op) :  mesh(m), ebhl(eb), turning_point(tp), obstacle_middle_point(op){ init();}
        EBHL_query_v2(EBHL_query_v2 const &) = delete;
        void operator=(EBHL_query_v2 const &x) = delete;

        ~ EBHL_query_v2()
        {
            shp_lab.Free();
        }
        void init(){
            timer = warthog::timer();
            polygons_check_list =  std::vector<int> (mesh->mesh_polygons.size());
            start_visibility_times_flag = std::vector<int> (turning_point.size());
            goal_visibility_times_flag = std::vector<int> (turning_point.size());
            start_visibility_cache = std::vector<bool> (turning_point.size());
            goal_visibility_cache = std::vector<bool> (turning_point.size());
            visible_search_id = 0;
        }

        void load_hub_label_on_convext_vertex( const string& label_file, const string& order_file){
            //load hub label
            shp_lab.load_labels(label_file.c_str());
            ifstream order_ifs(order_file.c_str());
            shp_inv =  vector<int>(numOfVertices);
            shp_rank = vector<int>(numOfVertices);
            for(int i = 0; i < numOfVertices; ++i){
                int tv;
                order_ifs >> tv;
                shp_rank[tv] = i;
                shp_inv[i] = tv;
            }

        }



        bool checkNextPolygon(int poly_id_1, int poly_id_2,const polyanya::Point& start_p, const polyanya::Point& end_p, const polyanya::PointLocation& end_p_pl){
            if(poly_id_2 == end_p_pl.poly1){
                assert(poly_id_2>=0 && poly_id_2 < polygons_check_list.size());

                //polygons_check_list[poly_id_2] = visible_search_id;
                polygons_check_list[poly_id_2] = visible_search_id;
                return true;

            }
            assert(poly_id_2>=0 && poly_id_2 < polygons_check_list.size());
            if( polygons_check_list[poly_id_2] == visible_search_id){
                polygons_check_list[poly_id_2] = visible_search_id;
                return false;

            }

            polygons_check_list[poly_id_2] = visible_search_id;



            const polyanya::Polygon& polygon = mesh->mesh_polygons[poly_id_2];
            const std::vector<int>& P = polygon.polygons;
            const std::vector<int>& V = polygon.vertices;
            const int N = (int) V.size();
            int last_vertex = V.back();
            for (int i = 0; i < N; i++)
            {
                const int this_vertex = V[i];
                assert(i < P.size());
                if (P[i] == poly_id_1)
                {
                    // The interval we're going to generate is the same as our

                    // current one, so skip it.

                    last_vertex = this_vertex;

                    continue;

                }
                assert(this_vertex >= 0 && this_vertex < mesh->mesh_vertices.size());
                assert(last_vertex >= 0 && last_vertex < mesh->mesh_vertices.size());
                const polyanya::Point& left = mesh->mesh_vertices[this_vertex].p,
                        right = mesh->mesh_vertices[last_vertex].p;
                if(doIntersect(left, right, start_p, end_p))
                {
                    if(left == end_p || right ==end_p){
                        return true;
                    }else  if(P[i] != -1 ){
                        if(checkNextPolygon(poly_id_2, P[i] , start_p, end_p, end_p_pl)){
                            return true;
                        }
                    }
                }
                last_vertex = this_vertex;
            }
            return false;
        }



        bool isVisible(const polyanya::Point& start_p, const polyanya::Point& end_p, const polyanya::PointLocation& start_p_pl, const polyanya::PointLocation& end_p_pl){
            visible_search_id ++;
            return checkNextPolygon(-1, start_p_pl.poly1, start_p, end_p, end_p_pl);

        }




        bool is_start_or_target_ambiguous(){
            const PointLocation& start_pl = get_point_location_in_search(start, mesh, false);
            const PointLocation& goal_pl = get_point_location_in_search(goal, mesh, false);
            if(start_pl.type == PointLocation::ON_CORNER_VERTEX_AMBIG){
                return true;
            }
            if(goal_pl.type == PointLocation::ON_CORNER_VERTEX_AMBIG){
                return true;
            }
            return false;
        }

        void set_start_goal(Point s, Point g)
        {
            start = s;
            goal = g;
        }

        //finds the shortest distance for a given start and target
        bool search(){
            start_index  = ebhl->get_grid_index(start);
            goal_index  = ebhl->get_grid_index(goal);
            const PointLocation& start_pl = get_point_location_in_search(start, mesh, false);
            const PointLocation& goal_pl = get_point_location_in_search(goal, mesh, false);

            const unsigned& start_begin = ebhl->get_grid_labels_begin(start_index);
            const unsigned& target_begin = ebhl->get_grid_labels_begin(goal_index);
            distance = INF_WEIGHT;

            if(isVisible(start,goal,start_pl,goal_pl)){
                distance = start.distance(goal);
                return true;
            }
            //go through the hubs and labels stored in both start and target grid
            for (unsigned i = start_begin, j = target_begin; ; ) {
                int v1 =  ebhl->get_hub_id(i), v2 = ebhl->get_hub_id(j);
                if (v1 == numOfVertices) break;  // Sentinel
                //if v1 and v2 are the same hub
                if (v1 == v2) {
                    double lower_bound = ebhl->get_lower_bound(i) + ebhl->get_lower_bound(j);
                     if(lower_bound > distance ){
                         i += 2;
                         j += 2;
                         continue;
                     }
                    double start_distance = INF_WEIGHT;
                    //check through all the distances of stored convex vertex to hub
                    for(unsigned k = ebhl->get_convex_begin(i); k != ebhl->get_convex_end(i);++k ){
                        const Convex_vertices_label& start_cv = ebhl->get_convex_label(k);
                        double tmp_distance = start.distance(turning_point[start_cv.convex_vertex])
                                              + start_cv.distance;
                        //if current label is visible
                        if (start_cv.visibility) {
                            if (tmp_distance < start_distance) {
                                start_distance = tmp_distance;
                            }
                        } else {
                            //if not visible but a potential shortest distance
                            if (tmp_distance < start_distance) {
                                //checks whether visibility has been checked previously
                                if (start_visibility_times_flag[start_cv.convex_vertex] == visible_search_id) {
                                    if (start_visibility_cache[start_cv.convex_vertex]) {
                                        start_distance = tmp_distance;
                                    }
                                } else {
                                    //checks the visibility
                                    start_visibility_times_flag[start_cv.convex_vertex] = visible_search_id;
                                    start_visibility_cache[start_cv.convex_vertex] = ebhl->check_visiblity(start,
                                                                                                           obstacle_middle_point[start_cv.convex_vertex],
                                                                                                           turning_point[start_cv.convex_vertex],
                                                                                                           start_cv.convex_vertex
                                    );
                                    if (start_visibility_cache[start_cv.convex_vertex]) {
                                        start_distance = tmp_distance;
                                    }
                                }

                            }
                        }
                    }
                    double target_distance = INF_WEIGHT;
                    for(unsigned k = ebhl->get_convex_begin(j); k != ebhl->get_convex_end(j);++k ){
                        const Convex_vertices_label& target_cv = ebhl->get_convex_label(k);
                        double tmp_distance = goal.distance(turning_point[target_cv.convex_vertex])
                                              + target_cv.distance;
                        if (target_cv.visibility) {
                            if (tmp_distance < target_distance) {
                                target_distance = tmp_distance;
                            }
                        } else {
                            if (tmp_distance < target_distance) {
                                if (goal_visibility_times_flag[target_cv.convex_vertex] == visible_search_id) {
                                    if (goal_visibility_cache[target_cv.convex_vertex]) {
                                        target_distance = tmp_distance;
                                    }
                                } else {
                                    goal_visibility_times_flag[target_cv.convex_vertex] = visible_search_id;
                                    goal_visibility_cache[target_cv.convex_vertex] = ebhl->check_visiblity(goal,
                                                                                                           obstacle_middle_point[target_cv.convex_vertex],
                                                                                                           turning_point[target_cv.convex_vertex],
                                                                                                           target_cv.convex_vertex
                                    );
                                    if (goal_visibility_cache[target_cv.convex_vertex]) {
                                        target_distance = tmp_distance;
                                    }
                                }

                            }

                        }

                    }
                    double td = start_distance + target_distance;
                    if (td < distance) {
                        distance = td;
                    }
                    i += 2;
                    j += 2;
                }
                else {
                    i += v1 < v2 ? 2 : 0;
                    j += v1 > v2 ? 2 : 0;
                }
            }
            //if start and target is not reachable
            if(distance == INF_WEIGHT){
                distance = INF;
            }
            timer.stop();
            return true;

        }


        //finds the shortest distance along with the path for a given start and target
        bool search_path(){
            start_index  = ebhl->get_grid_index(start);
            goal_index  = ebhl->get_grid_index(goal);
            const PointLocation& start_pl = get_point_location_in_search(start, mesh, false);
            const PointLocation& goal_pl = get_point_location_in_search(goal, mesh, false);
            distance = INF_WEIGHT;
            best_hub_id = -1;
            if(isVisible(start,goal,start_pl,goal_pl)){
                distance = start.distance(goal);
                return true;
            }
            const unsigned& start_begin = ebhl->get_grid_labels_begin(start_index);
            const unsigned& target_begin = ebhl->get_grid_labels_begin(goal_index);
            //goes through all the labels stored in both start and target grids
            for (unsigned i = start_begin, j = target_begin; ; ) {
                int v1 =  ebhl->get_hub_id(i), v2 = ebhl->get_hub_id(j);
                if (v1 == numOfVertices) break;  // Sentinel
                //if v1 and v2 are the same hub
                if (v1 == v2) {
                    double lower_bound = ebhl->get_lower_bound(i) + ebhl->get_lower_bound(j);
                    if(lower_bound > distance ){
                        i += 2;
                        j += 2;
                        continue;
                    }
                    double start_distance = INF_WEIGHT;
                    int best_start_predecessor = -1;
                    int best_start_convex = -1;
                    //check through all the distances of stored convex vertex to hub
                    for(unsigned k = ebhl->get_convex_begin(i); k != ebhl->get_convex_end(i);++k ){
                        const Convex_vertices_label& start_cv = ebhl->get_convex_label(k);
                        double tmp_distance = start.distance(turning_point[start_cv.convex_vertex])
                                              + start_cv.distance;
                        //if the visibility is true
                        if (start_cv.visibility) {
                            if (tmp_distance < start_distance) {
                                start_distance = tmp_distance;
                                best_start_predecessor = start_cv.predecessor;
                                best_start_convex  = start_cv.convex_vertex;
                            }

                        } else {
                            //not visible but a potential shortest distance
                            if (tmp_distance < start_distance) {
                                //visibility has been checked previously
                                if (start_visibility_times_flag[start_cv.convex_vertex] == visible_search_id) {
                                    if (start_visibility_cache[start_cv.convex_vertex]) {
                                        start_distance = tmp_distance;
                                        best_start_predecessor = start_cv.predecessor;
                                        best_start_convex  = start_cv.convex_vertex;
                                    }
                                } else {
                                    //checks visibility on the fly
                                    start_visibility_times_flag[start_cv.convex_vertex] = visible_search_id;
                                    start_visibility_cache[start_cv.convex_vertex] = ebhl->check_visiblity(start,
                                                                                                           obstacle_middle_point[start_cv.convex_vertex],
                                                                                                           turning_point[start_cv.convex_vertex],
                                                                                                           start_cv.convex_vertex
                                    );

                                    if (start_visibility_cache[start_cv.convex_vertex]) {
                                        start_distance = tmp_distance;
                                        best_start_predecessor = start_cv.predecessor;
                                        best_start_convex  = start_cv.convex_vertex;
                                    }
                                }

                            }
                        }
                    }
                    double target_distance = INF_WEIGHT;
                    int best_target_predecessor = -1;
                    int best_target_convex = -1;
                    for(unsigned k = ebhl->get_convex_begin(j); k != ebhl->get_convex_end(j);++k ){
                        const Convex_vertices_label& target_cv = ebhl->get_convex_label(k);
                        double tmp_distance = goal.distance(turning_point[target_cv.convex_vertex])
                                              + target_cv.distance;

                        if (target_cv.visibility) {
                            if (tmp_distance < target_distance) {
                                target_distance = tmp_distance;
                                best_target_predecessor = target_cv.predecessor;
                                best_target_convex = target_cv.convex_vertex;
                            }
                        } else {
                            if (tmp_distance < target_distance) {
                                if (goal_visibility_times_flag[target_cv.convex_vertex] == visible_search_id) {
                                    if (goal_visibility_cache[target_cv.convex_vertex]) {
                                        target_distance = tmp_distance;
                                        best_target_predecessor = target_cv.predecessor;
                                        best_target_convex = target_cv.convex_vertex;
                                    }
                                } else {
                                    goal_visibility_times_flag[target_cv.convex_vertex] = visible_search_id;
                                    goal_visibility_cache[target_cv.convex_vertex] = ebhl->check_visiblity(goal,
                                                                                                           obstacle_middle_point[target_cv.convex_vertex],
                                                                                                           turning_point[target_cv.convex_vertex],
                                                                                                           target_cv.convex_vertex
                                    );

                                    if (goal_visibility_cache[target_cv.convex_vertex]) {
                                        target_distance = tmp_distance;
                                        best_target_predecessor = target_cv.predecessor;
                                        best_target_convex = target_cv.convex_vertex;
                                    }
                                }

                            }

                        }

                    }
                    double td = start_distance + target_distance;
                    if (td < distance) {
                        distance = td;
                        best_hub_id = v1;
                        target_predecessor = best_target_predecessor;
                        target_convex = best_target_convex;
                        start_predecessor = best_start_predecessor;
                        start_convex = best_start_convex;
                    }
                    i += 2;
                    j += 2;
                }
                else {
                    i += v1 < v2 ? 2 : 0;
                    j += v1 > v2 ? 2 : 0;
                }
            }
            //start and target is not reachable
            if(distance == INF_WEIGHT){
                distance = INF;
            }
            timer.stop();
            return true;
        }

        double get_cost(){
            return distance;
        }

        double get_search_micro()
        {
            return timer.elapsed_time_micro();
        }

        double get_search_nano()
        {
            return timer.elapsed_time_nano();
        }

        //retrieves path for given start and target if reachable
        vector<Point> get_path(){
            shortest_path.clear();
            if(distance != -1){
                if(best_hub_id == -1){
                    shortest_path.push_back(start);
                    shortest_path.push_back(goal);
                }else{
                    shortest_path.push_back(start);
                    shp_lab.query_path(start_convex,start_predecessor,target_convex,target_predecessor,best_hub_id,
                                       shp_rank, shp_inv, turning_point,shortest_path);
                    shortest_path.push_back(goal);
                }

                Point current = shortest_path[0];
                double test_distance = 0;
                for (int i = 1; i < shortest_path.size(); i++) {
                    test_distance = test_distance  + shortest_path[i].distance(current);
                    current =  shortest_path[i];
                }
                if (fabs(test_distance - distance)>EPSILON) {
                    std::cout << distance << std::endl;
                    std::cout << test_distance << std::endl;
                }

            }
            return shortest_path;
        }

    };


}

#endif //LEARNED_HEURISTIC_EBHL_QUERY_V2_H

