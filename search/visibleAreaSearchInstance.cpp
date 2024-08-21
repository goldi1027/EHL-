/*
 * Modified version of search instance class from polyanya
 * Changes the search of visible nodes into visible space
 */

#include "visibleAreaSearchInstance.h"
#include "expansion.h"
#include "geometry.h"
#include "searchnode.h"
#include "successor.h"
#include "vertex.h"
#include "mesh.h"
#include "point.h"
#include "consts.h"
#include <queue>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <unordered_set>
namespace bg = boost::geometry;
using namespace std;
namespace polyanya {

    int visibleAreaSearchInstance::succ_to_node(
            SearchNodePtr parent, Successor *successors, int num_succ,
            SearchNode *nodes
    ) {
        assert(mesh != nullptr);
#ifndef NDEBUG
        if (verbose) {
            for (int i = 0; i < num_succ; i++) {
                const Successor &succ = successors[i];
                std::cerr << "received successors: ";
                std::cerr << succ;
                std::cerr << std::endl;
            }
        }
#endif

        const Polygon &polygon = mesh->mesh_polygons[parent->next_polygon];
        const std::vector<int> &V = polygon.vertices;
        const std::vector<int> &P = polygon.polygons;

        double right_g = -1, left_g = -1;

        int out = 0;

        for (int i = 0; i < num_succ; i++) {
            const Successor &succ = successors[i];
            const int next_polygon = P[succ.poly_left_ind];
            //turn off here to generate vertices on obstacle corner

            //if (next_polygon == -1)
            //{
            //    continue;
            //}

            // If the successor we're about to push pushes into a one-way polygon,
            // and the polygon isn't the end polygon, just continue.

            //turn this off if we want to find all vetices.
//            if (mesh->mesh_polygons[next_polygon].is_one_way &&
//                next_polygon != end_polygon)
//            {
//                continue;
//            }
            const int left_vertex = V[succ.poly_left_ind];
            const int right_vertex = succ.poly_left_ind ?
                                     V[succ.poly_left_ind - 1] :
                                     V.back();

            // Note that g is evaluated twice here. (But this is a lambda!)
            // Always try to precompute before using this macro.
            // We implicitly set h to be zero and let search() update it.
            const auto p = [&](const int root, const double g) {
                if (root != -1) {
                    assert(root >= 0 && root < (int) root_g_values.size());
                    // Can POSSIBLY prune?
                    if (root_search_ids[root] != search_id) {
                        // First time reaching root
                        root_search_ids[root] = search_id;
                        root_g_values[root] = g;
                    } else {
                        // We've been here before!
                        // Check whether we've done better.
                        if (root_g_values[root] + EPSILON < g) {
                            // We've done better!
                            return;
                        } else {
                            // This is better.
                            root_g_values[root] = g;
                        }

                    }

                }
                nodes[out++] = {nullptr, root, succ.left, succ.right, left_vertex,
                                right_vertex, next_polygon, g, g};
            };

            const Point &parent_root = (parent->root == -1 ?
                                        start :
                                        mesh->mesh_vertices[parent->root].p);
#define get_g(new_root) parent->g + parent_root.distance(new_root)

            switch (succ.type) {
                case Successor::RIGHT_NON_OBSERVABLE:
                    if (right_g == -1) {
                        right_g = get_g(parent->right);
                    }
                    p(parent->right_vertex, right_g);
                    break;

                case Successor::OBSERVABLE:
                    p(parent->root, parent->g);
                    break;

                case Successor::LEFT_NON_OBSERVABLE:
                    if (left_g == -1) {
                        left_g = get_g(parent->left);
                    }
                    p(parent->left_vertex, left_g);
                    break;

                default:
                    assert(false);
                    break;
            }
#undef get_h
#undef get_g
        }

        return out;
    }

#define root_to_point(root) ((root) == -1 ? start : mesh->mesh_vertices[root].p)




    double get_angles_between_points(const Point& p1, const Point& p2, const Point& p3) {
        double AB = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
        double BC = sqrt(pow(p2.x - p3.x, 2) + pow(p2.y - p3.y, 2));
        double AC = sqrt(pow(p3.x - p1.x, 2) + pow(p3.y - p1.y, 2));
        double value = (BC * BC + AB * AB - AC * AC) / (2 * BC * AB);
        value = value < -1 ? -1 : value;
        double raidus = acos(value);
        double degree = raidus * 180 / PI;
        if (get_orientation(p1, p2, p3) == Orientation::CCW) {
            degree = 360 - degree;
        }
        return degree;
    }

//
//    double get_angles_between_points(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y) {
//
//        double AB = sqrt(pow(p2x - p1x, 2) + pow(p2y - p1y, 2));
//        double BC = sqrt(pow(p2x - p3x, 2) + pow(p2y - p3y, 2));
//        double AC = sqrt(pow(p3x - p1x, 2) + pow(p3y - p1y, 2));
//        double raidus = acos((BC * BC + AB * AB - AC * AC) / (2 * BC * AB));
//        double degree = raidus * 180 / PI;
//        auto a = get_orientation(Point{p1x, p1y}, Point{p2x, p2y}, Point{p3x, p3y});
//        if (get_orientation(Point{p1x, p1y}, Point{p2x, p2y}, Point{p3x, p3y}) == Orientation::CCW) {
//            degree = 360 - degree;
//        }
//        return degree;
//    }



    void visibleAreaSearchInstance::gen_initial_nodes_for_search_visibility_area(const PointLocation& start_pl)
    {
        // {parent, root, left, right, next_polygon, right_vertex, f, g}
        // be VERY lazy and abuse how our function expands collinear search nodes
        // if right_vertex is not valid, it will generate EVERYTHING
        // and we can set right_vertex if we want to omit generating an interval.
        const PointLocation pl = start_pl;

#define get_lazy(next, left, right) new (node_pool->allocate()) SearchNode \
        {nullptr, -1, start, start, left, right, next, 0, 0}

#define v(vertex) mesh->mesh_vertices[vertex]

        const auto push_lazy = [&](SearchNodePtr lazy)
        {
            const int poly = lazy->next_polygon;
            if (poly == -1)
            {
                return;
            }
            // iterate over poly, throwing away vertices if needed
            const std::vector<int>& vertices =
                    mesh->mesh_polygons[poly].vertices;
            Successor* successors = new Successor [vertices.size()];
            int last_vertex = vertices.back();
            int num_succ = 0;
            for (int i = 0; i < (int) vertices.size(); i++)
            {
                const int vertex = vertices[i];
                const Vertex& v = mesh->mesh_vertices[vertex];
                if (vertex == lazy->right_vertex ||
                    last_vertex == lazy->left_vertex)
                {
                    last_vertex = vertex;
                    continue;
                }
                successors[num_succ++] =
                        {Successor::OBSERVABLE, v(vertex).p,
                         v(last_vertex).p, i};
                last_vertex = vertex;
            }
            SearchNode* nodes = new SearchNode [num_succ];
            const int num_nodes = succ_to_node(lazy, successors,
                                               num_succ, nodes);
            delete[] successors;
            for (int i = 0; i < num_nodes; i++)
            {
                SearchNodePtr n = new (node_pool->allocate())
                        SearchNode(nodes[i]);
                const Point& n_root = (n->root == -1 ? start :
                                       mesh->mesh_vertices[n->root].p);
                n->f += get_interval_heuristic(n_root, n->left, n->right);
                n->parent = lazy;
                open_list.push(n);
            }
            delete[] nodes;
            nodes_generated += num_nodes;
            nodes_pushed += num_nodes;
        };
        //since all the node on conrner_vertex do not handle other type
        SearchNodePtr lazy = get_lazy(pl.poly1, -1, -1);
        push_lazy(lazy);
        nodes_generated++;
#undef v
#undef get_lazy
    }

    void visibleAreaSearchInstance::gen_initial_nodes_for_search_non_taut_visibility_area(const PointLocation& start_pl)
    {
        // {parent, root, left, right, next_polygon, right_vertex, f, g}
        // be VERY lazy and abuse how our function expands collinear search nodes
        // if right_vertex is not valid, it will generate EVERYTHING
        // and we can set right_vertex if we want to omit generating an interval.
        const PointLocation pl = start_pl;
        Point obstacle_point1 =mesh->mesh_vertices[ mesh->mesh_vertices[start_id].obstacle_edge[0]].p;
        Point obstacle_point2 =mesh->mesh_vertices[ mesh->mesh_vertices[start_id].obstacle_edge[1]].p;
        // always calcuate angle based on obstacle_point -> s -> other_vertices.
        Point first_obstacle_point = get_orientation(obstacle_point1,start,obstacle_point2) == Orientation::CCW ? obstacle_point1 : obstacle_point2;
        Point second_obstacle_point = obstacle_point1 == first_obstacle_point ?  obstacle_point2 : obstacle_point1;
        //add small epsilon
        double min_angle = get_angles_between_points(first_obstacle_point, start, second_obstacle_point) - 180  - EPSILON;
        double max_angle = 180 + EPSILON;




#define get_lazy(next, left, right) new (node_pool->allocate()) SearchNode \
        {nullptr, -1, start, start, left, right, next, 0, 0}

#define v(vertex) mesh->mesh_vertices[vertex]

        const auto push_lazy = [&](SearchNodePtr lazy)
        {
            const int poly = lazy->next_polygon;
            if (poly == -1)
            {
                return;
            }
            // iterate over poly, throwing away vertices if needed
            const std::vector<int>& vertices =
                    mesh->mesh_polygons[poly].vertices;
            Successor* successors = new Successor [vertices.size()];
            int last_vertex = vertices.back();
            int num_succ = 0;
            for (int i = 0; i < (int) vertices.size(); i++)
            {
                const int vertex = vertices[i];
                const Vertex& v = mesh->mesh_vertices[vertex];
                if (vertex == lazy->right_vertex ||
                    last_vertex == lazy->left_vertex)
                {
                    last_vertex = vertex;
                    continue;
                }
                successors[num_succ++] =
                        {Successor::OBSERVABLE, v(vertex).p,
                         v(last_vertex).p, i};
                last_vertex = vertex;
            }
            SearchNode* nodes = new SearchNode [num_succ];
            const int num_nodes = succ_to_node(lazy, successors,
                                               num_succ, nodes);
            delete[] successors;
            for (int i = 0; i < num_nodes; i++)
            {
                SearchNodePtr n = new (node_pool->allocate())
                        SearchNode(nodes[i]);
                const Point& n_root = (n->root == -1 ? start :
                                       mesh->mesh_vertices[n->root].p);
                if(is_collinear(start, n->left,n->right )){
                    continue;
                }

                double left_angle = get_angles_between_points(first_obstacle_point,start,n->left);
                double right_angle = get_angles_between_points(first_obstacle_point,start,n->right);

                if(min_angle <= left_angle && left_angle<= max_angle && min_angle <= right_angle && right_angle<= max_angle){
                    // successor is compeleted within non-taut region
                    continue;
                }else if((min_angle <= left_angle && left_angle<= max_angle) && !(min_angle <= right_angle && right_angle<= max_angle) ){
                    //left vertex within the non-taut region, right is outisde the non-taut region
                    Point intersection = line_intersect(second_obstacle_point,start, n->left,n->right );
                    n->left = intersection;
                    n->f += get_interval_heuristic(n_root, n->left, n->right);
                    //std::cout<<"Generate successors: " << n->left << n->right;
                    n->parent = lazy;
                    open_list.push(n);
                }else if (!(min_angle <= left_angle && left_angle<= max_angle) && (min_angle <= right_angle && right_angle<= max_angle) ){
                    //right vertex within the non-taut region, left is outisde the non-taut region
                    Point intersection = line_intersect(first_obstacle_point,start, n->left,n->right );
                    n->right = intersection;
                    n->f += get_interval_heuristic(n_root, n->left, n->right);
                    //std::cout<<"Generate successors: " << n->left << n->right;
                    n->parent = lazy;
                    open_list.push(n);
                }else{
                    // use mid point to test whether a successors is intersect with non-taut region;
                    Point mid_point = {(first_obstacle_point.x+ second_obstacle_point.x)/2,
                                       (first_obstacle_point.y+ second_obstacle_point.y)/2};

                    if(raytTolineIntersection(mid_point,start, n->left,n->right)){
                        // this part may contain error, check this with biggger map.
                        // generate left -> left_intersection;
                        Point left_intersection = line_intersect(first_obstacle_point,start, n->left,n->right );
                        n->right = left_intersection;
                        n->f += get_interval_heuristic(n_root, n->left, n->right);
                        n->parent = lazy;
                        open_list.push(n);
                        //std::cout<<"Generate successors: " << n->left << n->right;
                        // generate right_intersection -> right
                        SearchNodePtr n2 = new (node_pool->allocate())
                                SearchNode(nodes[i]);
                        Point right_intersection = line_intersect(second_obstacle_point,start, n->left,n->right );
                        n2->left = right_intersection;
                        n2->f += get_interval_heuristic(n_root, n2->left, n2->right);
                        n2->parent = lazy;
                        open_list.push(n2);
                        //std::cout<<"Generate successors: " << n2->left << n2->right;

                    }else{
                        n->f += get_interval_heuristic(n_root, n->left, n->right);
                        n->parent = lazy;
                        open_list.push(n);
                        //std::cout<<"Generate successors: " << n->left << n->right;
                    }
                }
            }
            delete[] nodes;
            nodes_generated += num_nodes;
            nodes_pushed += num_nodes;
        };
        //since all the node on conrner_vertex do not handle other type
        //initalize all adjacent polygons.
        for(  int& p_id :mesh->mesh_vertices[start_id].polygons){
            if(p_id == -1 ){
                continue;
            }
            SearchNodePtr lazy = get_lazy(p_id, -1, -1);
            push_lazy(lazy);
            nodes_generated++;
        }

#undef v
#undef get_lazy
    }




    /*
     * Search for visible area
     * Pruning 1 is applied which splits the visible area into two non-taut regions
     */
    bool visibleAreaSearchInstance::search_non_taut_visible_area(int s,PointLocation start_pl,polygon_type& poly1, polygon_type& poly2)
    {

        timer.start();
        vector<tuple<Point,Point>> successors_on_obstacles;
        start = mesh->mesh_vertices[s].p;
        start_id = s;
        init_search_non_taut_visible_area(start_pl);
        if (mesh == nullptr)
        {
            timer.stop();
            return false;
        }

        while (!open_list.empty())
        {

            SearchNodePtr node = open_list.top(); open_list.pop();
            if(node->next_polygon == -1){
                if(get_orientation(start,node->left,node->right) != Orientation::COLLINEAR){
                    successors_on_obstacles.emplace_back(node->left,node->right);
                }
                continue;
            }
            nodes_popped++;
            // We will never update our root list here.
            const int root = node->root;
            if(root != -1 && root != s){
                std::cout<<"root changed"<<std::endl;
            }
            int num_nodes = 1;
            search_nodes_to_push[0] = *node;

            SearchNode cur_node = search_nodes_to_push[0];
            int num_succ = get_observable_successors(cur_node, start, *mesh,
                                                     search_successors);
            successor_calls++;
            num_nodes = succ_to_node(&cur_node, search_successors,
                                     num_succ, search_nodes_to_push);

            for (int i = 0; i < num_nodes; i++)
            {
                // We need to update the h value before we push!
                const SearchNodePtr n = new (node_pool->allocate())
                        SearchNode(search_nodes_to_push[i]);
                const Point& n_root = (n->root == -1 ? start :
                                       mesh->mesh_vertices[n->root].p);
                n->f += get_interval_heuristic(n_root, n->left, n->right);


                // This node's parent should be nullptr, so we should set it.
                n->parent = node;
                open_list.push(n);
            }
            nodes_generated += num_nodes;
            nodes_pushed += num_nodes;
        }

        Point obstacle_point1 =mesh->mesh_vertices[ mesh->mesh_vertices[s].obstacle_edge[0]].p;
        Point obstacle_point2 =mesh->mesh_vertices[ mesh->mesh_vertices[s].obstacle_edge[1]].p;
        Point first_obstacle_point = get_orientation(obstacle_point1,start,obstacle_point2) == Orientation::CCW ? obstacle_point1 : obstacle_point2;
        Point second_obstacle_point = obstacle_point1 == first_obstacle_point?  obstacle_point2 : obstacle_point1;

        double min_angle = get_angles_between_points(first_obstacle_point, start, second_obstacle_point) - 180 ;
        double max_angle = 180;


        Point start_point = start;
        std::sort(std::begin(successors_on_obstacles),
                  std::end(successors_on_obstacles),
                  [ & first_obstacle_point, & start_point]( const tuple<Point,Point>&  lvalue, const tuple<Point,Point>& rvalue) {
                      double angle_1 =  get_angles_between_points(first_obstacle_point, start_point,get<1>(lvalue));
                      double angle_2 =  get_angles_between_points(first_obstacle_point, start_point,get<1>(rvalue));
                      return angle_1 < angle_2; });

        bg::append( poly1.outer(), point_type{start.x,start.y});
        unsigned second_polygon_start_index = 0 ;
        for(unsigned i = 0 ;  i < successors_on_obstacles.size(); i ++){
            const auto& successor = successors_on_obstacles[i];
            if(get_angles_between_points(first_obstacle_point, start_point, get<1>(successor)) > min_angle){
                second_polygon_start_index = i ;
                break;
            }
            if( !(poly1.outer().back().x() == get<1>(successor).x && poly1.outer().back().y() == get<1>(successor).y)){
                bg::append( poly1.outer(), point_type{get<1>(successor).x,get<1>(successor).y});
            }
            bg::append( poly1.outer(), point_type{get<0>(successor).x,get<0>(successor).y});
        }
        bg::append( poly1.outer(), point_type{start.x,start.y});
        // sort by CCW;
        bg::correct(poly1);

        if(poly1.outer().size()>3){
            polygon_type poly_refined1;
            bg::append( poly_refined1.outer(), point_type{start.x,start.y});
            for(unsigned i  = 1 ; i < poly1.outer().size()-1; i++){

                if(!is_collinear(Point{poly_refined1.outer().back().x(),poly_refined1.outer().back().y()},
                                 Point{poly1.outer()[i].x(),poly1.outer()[i].y()},
                                 Point{poly1.outer()[i+1].x(),poly1.outer()[i+1].y()})){
                    bg::append( poly_refined1.outer(), poly1.outer()[i]);
                }
            }
            bg::append( poly_refined1.outer(), point_type{start.x,start.y});
            poly1 = poly_refined1;
        }

        if(second_polygon_start_index == 0 || second_polygon_start_index == successors_on_obstacles.size()-1 ){
            std::cout<<"only find one polygon"<< s << std::endl;
        }

        bg::append( poly2.outer(), point_type{start.x,start.y});
        for(unsigned i = second_polygon_start_index ;  i <  successors_on_obstacles.size(); i ++){
            const auto& successor = successors_on_obstacles[i];
            if( !(poly2.outer().back().x() == get<1>(successor).x && poly2.outer().back().y() == get<1>(successor).y)){
                bg::append( poly2.outer(), point_type{get<1>(successor).x,get<1>(successor).y});
            }
            bg::append( poly2.outer(), point_type{get<0>(successor).x,get<0>(successor).y});
        }
        // sort by CCW;
        bg::correct(poly2);

        if(poly2.outer().size()>3){
            polygon_type poly_refined2;
            bg::append( poly_refined2.outer(), point_type{start.x,start.y});
            for(unsigned i  = 1 ; i < poly2.outer().size()-1; i++){

                if(!is_collinear(Point{poly_refined2.outer().back().x(),poly_refined2.outer().back().y()},
                                 Point{poly2.outer()[i].x(),poly2.outer()[i].y()},
                                 Point{poly2.outer()[i+1].x(),poly2.outer()[i+1].y()})){
                    bg::append( poly_refined2.outer(), poly2.outer()[i]);
                }
            }
            bg::append( poly_refined2.outer(), point_type{start.x,start.y});
            poly2 = poly_refined2;
        }
//
//        polygon_type poly_refined2;
//        bg::append( poly_refined2.outer(), point_type{start.x,start.y});
//        for(unsigned i  = 1 ; i < poly1.outer().size()-1; i++){
//
//            if(!is_collinear(Point{poly_refined2.outer().back().x(),poly_refined2.outer().back().y()},
//                             Point{poly2.outer()[i].x(),poly2.outer()[i].y()},
//                             Point{poly2.outer()[i+1].x(),poly2.outer()[i+1].y()})){
//                bg::append( poly_refined2.outer(), poly2.outer()[i]);
//            }
//        }
//        bg::append( poly_refined2.outer(), point_type{start.x,start.y});
//        bg::correct(poly_refined2);
//        poly2 = poly_refined2;

//        for (const auto& p : polygon1 ){
//            std::cout<< p << std::endl;
//            std::cout<< get_angles_between_points(first_obstacle_point, start, p) << std::endl;
//        }
//
//        for (const auto& p : polygon2 ){
//            std::cout<< p << std::endl;
//            std::cout<< get_angles_between_points(first_obstacle_point, start, p) << std::endl;
//        }
        timer.stop();
        return false;
    }


    bool visibleAreaSearchInstance::search_visible_area(int s,PointLocation start_pl, polygon_type& poly)
    {
        timer.start();
        vector<tuple<Point,Point>> successors_on_obstacles;
        start = mesh->mesh_vertices[s].p;
        init_search_visible_area(start_pl);
        if (mesh == nullptr)
        {
            timer.stop();
            return false;
        }

        while (!open_list.empty())
        {

            SearchNodePtr node = open_list.top(); open_list.pop();
            if(node->next_polygon == -1){
                if(get_orientation(start,node->left,node->right) != Orientation::COLLINEAR){
                    successors_on_obstacles.emplace_back(node->left,node->right);
                }
                continue;
            }
            nodes_popped++;
            // We will never update our root list here.
            const int root = node->root;
            if(root != -1 && root != s){
                std::cout<<"root changed"<<std::endl;
            }
            int num_nodes = 1;
            search_nodes_to_push[0] = *node;

            SearchNode cur_node = search_nodes_to_push[0];
            int num_succ = get_observable_successors(cur_node, start, *mesh,
                                                     search_successors);
            successor_calls++;
            num_nodes = succ_to_node(&cur_node, search_successors,
                                     num_succ, search_nodes_to_push);

            for (int i = 0; i < num_nodes; i++)
            {
                // We need to update the h value before we push!
                const SearchNodePtr n = new (node_pool->allocate())
                        SearchNode(search_nodes_to_push[i]);
                const Point& n_root = (n->root == -1 ? start :
                                       mesh->mesh_vertices[n->root].p);
                n->f += get_interval_heuristic(n_root, n->left, n->right);


                // This node's parent should be nullptr, so we should set it.
                n->parent = node;
                open_list.push(n);
            }
            nodes_generated += num_nodes;
            nodes_pushed += num_nodes;
        }

        Point obstacle_point1 =mesh->mesh_vertices[ mesh->mesh_vertices[s].obstacle_edge[0]].p;
        Point obstacle_point2 =mesh->mesh_vertices[ mesh->mesh_vertices[s].obstacle_edge[1]].p;
        Point obstacle_point = get_orientation(obstacle_point1,start,obstacle_point2) == Orientation::CCW ? obstacle_point1 : obstacle_point2;

        Point start_point = start;
        std::sort(std::begin(successors_on_obstacles),
                  std::end(successors_on_obstacles),
                  [ & obstacle_point, & start_point]( const tuple<Point,Point>&  lvalue, const tuple<Point,Point>& rvalue) {
                      double angle_1 =  get_angles_between_points(obstacle_point, start_point,get<1>(lvalue));
                      double angle_2 =  get_angles_between_points(obstacle_point, start_point,get<1>(rvalue));
                      return angle_1 < angle_2; });

        bg::append( poly.outer(), point_type{start.x,start.y});
        for(const auto& successor : successors_on_obstacles ){
            if( !(poly.outer().back().x() == get<1>(successor).x && poly.outer().back().y() == get<1>(successor).y)){
                bg::append( poly.outer(), point_type{get<1>(successor).x,get<1>(successor).y});
            }
            bg::append( poly.outer(), point_type{get<0>(successor).x,get<0>(successor).y});
        }
        bg::correct(poly);
//        for (const auto& p : polygon ){
//            std::cout<< p << std::endl;
//            std::cout<< get_angles_between_points(obstacle_point.x,obstacle_point.y, start_point.x, start_point.y,p.x,p.y) << std::endl;
//        }
        timer.stop();
        return false;
    }










    void visibleAreaSearchInstance::print_node(SearchNodePtr node, std::ostream& outfile)
    {
        outfile << "root=" << root_to_point(node->root) << "; left=" << node->left
                << "; right=" << node->right << "; f=" << node->f << ", g="
                << node->g;
        /*
        outfile << "; col=" << [&]() -> std::string
                {
                    switch (node->col_type)
                    {
                        case SearchNode::NOT:
                            return "NOT";
                        case SearchNode::RIGHT:
                            return "RIGHT";
                        case SearchNode::LEFT:
                            return "LEFT";
                        case SearchNode::LAZY:
                            return "LAZY";
                        default:
                            return "";
                    }
                }();
        */
    }




#undef root_to_point

}