/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */
#include "visibleSearchInstance.h"
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
using namespace std;
namespace polyanya {

    int visibleSearchInstance::succ_to_node(
            SearchNodePtr parent, Successor* successors, int num_succ,
            SearchNode* nodes
    )
    {
        assert(mesh != nullptr);
#ifndef NDEBUG
        if (verbose)
        {    for (int i = 0; i < num_succ; i++) {
                const Successor &succ = successors[i];
                std::cerr << "received successors: ";
                std::cerr << succ;
                std::cerr << std::endl;
            }
        }
#endif

        const Polygon& polygon = mesh->mesh_polygons[parent->next_polygon];
        const std::vector<int>& V = polygon.vertices;
        const std::vector<int>& P = polygon.polygons;

        double right_g = -1, left_g = -1;

        int out = 0;




        for (int i = 0; i < num_succ; i++)
        {

            const Successor& succ = successors[i];
            const int next_polygon = P[succ.poly_left_ind];
            if (next_polygon == -1)
            {
                continue;
            }


            // If the successor we're about to push pushes into a one-way polygon,
            // and the polygon isn't the end polygon, just continue.

            //turn this off if we want to find all vetices.
//            if (mesh->mesh_polygons[next_polygon].is_one_way &&
//                next_polygon != end_polygon)
//            {
//                continue;
//            }
            const int left_vertex  = V[succ.poly_left_ind];
            const int right_vertex = succ.poly_left_ind ?
                                     V[succ.poly_left_ind - 1] :
                                     V.back();

            // Note that g is evaluated twice here. (But this is a lambda!)
            // Always try to precompute before using this macro.
            // We implicitly set h to be zero and let search() update it.
            const auto p = [&](const int root, const double g)
            {
                if (root != -1)
                {
                    assert(root >= 0 && root < (int) root_g_values.size());
                    // Can POSSIBLY prune?
                    if (root_search_ids[root] != search_id)
                    {
                        // First time reaching root
                        root_search_ids[root] = search_id;
                        root_g_values[root] = g;
                    }
                    else
                    {
                        // We've been here before!
                        // Check whether we've done better.
                        if (root_g_values[root] + EPSILON <  g  )
                        {
                            // We've done better!
                            return;
                        }else{
                            // This is better.
                            root_g_values[root] = g;
                        }

                    }

                }
                nodes[out++] = {nullptr, root, succ.left, succ.right, left_vertex,
                                right_vertex, next_polygon, g, g};
            };

            const Point& parent_root = (parent->root == -1 ?
                                        start :
                                        mesh->mesh_vertices[parent->root].p);
#define get_g(new_root) parent->g + parent_root.distance(new_root)

            switch (succ.type)
            {
                case Successor::RIGHT_NON_OBSERVABLE:
                    if (right_g == -1)
                    {
                        right_g = get_g(parent->right);
                    }
                    p(parent->right_vertex, right_g);
                    break;

                case Successor::OBSERVABLE:
                    p(parent->root, parent->g);
                    break;

                case Successor::LEFT_NON_OBSERVABLE:
                    if (left_g == -1)
                    {
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


    void visibleSearchInstance::gen_initial_nodes_for_visiblity_search(const PointLocation& start_pl, vector<int>& visible_vertex)
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
                if(v.p != start){
                    if(visible_vertices_list[vertex] != vertice_search_id){
                        //not retrieved before;
                        if(!v.is_ambig && v.is_turning_vertex) {
                            visible_vertex.push_back(vertex);
//                            visible_vertex.insert(vertex);
                        }
                        visible_vertices_list[vertex] = vertice_search_id;
                        degree ++;
                    }

                }
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



    bool visibleSearchInstance::search_visible_vertices(int s,PointLocation start_pl,vector<int>& visible_vertex)
    {
        timer.start();
        start = mesh->mesh_vertices[s].p;
        init_search(start_pl,visible_vertex);

        if (mesh == nullptr)
        {
            timer.stop();
            return false;
        }

        while (!open_list.empty())
        {

            SearchNodePtr node = open_list.top(); open_list.pop();

            const Vertex&  left_v = mesh->mesh_vertices[ node->left_vertex];
            const Vertex&  right_v = mesh->mesh_vertices[ node->right_vertex];
            if(node->left == left_v.p  && left_v.p != start) {
                if (visible_vertices_list[node->left_vertex] != vertice_search_id) {
                    if (left_v.is_turning_vertex && !left_v.is_ambig ) {
                        visible_vertex.push_back(node->left_vertex);
                    }
                    degree ++;
                    visible_vertices_list[node->left_vertex] = vertice_search_id;
                }
            }
            if(node->right == right_v.p && right_v.p != start){
                if(visible_vertices_list[node->right_vertex] != vertice_search_id) {
                    if (right_v.is_turning_vertex&&!right_v.is_ambig) {
                        visible_vertex.push_back(node->right_vertex);
                    }
                    degree ++;
                    visible_vertices_list[node->right_vertex] = vertice_search_id;
                }
            }

            nodes_popped++;
            const int next_poly = node->next_polygon;
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
//                if(get_orientation(start,search_nodes_to_push[i].left,search_nodes_to_push[i].right) == Orientation::COLLINEAR){
//                    if(start != search_nodes_to_push[i].left && search_nodes_to_push[i].right !=start){
//                        continue;
//                    }
//                }
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

        timer.stop();
        return false;
    }



    void visibleSearchInstance::print_node(SearchNodePtr node, std::ostream& outfile)
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