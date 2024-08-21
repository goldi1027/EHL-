/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */
#pragma once
#include "searchnode.h"
#include "successor.h"
#include "mesh.h"
#include "point.h"
#include "cpool.h"
#include "timer.h"
#include <queue>
#include <vector>
#include <ctime>
#include <chrono>
#include <algorithm>


namespace polyanya
{

    template<typename T, typename Compare = std::greater<T> >
    struct PointerComp
    {
        bool operator()(const T* x,
                        const T* y) const
        {
            return Compare()(*x, *y);
        }
    };

    typedef Mesh* MeshPtr;
// Polyanya instance for point to point search
    class SearchInstance
    {
        typedef std::priority_queue<SearchNodePtr, std::vector<SearchNodePtr>,
                PointerComp<SearchNode> > pq;
    private:
        warthog::mem::cpool* node_pool;
        MeshPtr mesh;
        Point start, goal;

        SearchNodePtr final_node;
        int end_polygon; // set by init_search
        pq open_list;

        // Best g value for a specific vertex.
        std::vector<double> root_g_values;
        // Contains the current search id if the root has been reached by
        // the search.
        std::vector<int> root_search_ids;  // also used for root-level pruning

        int search_id;

        warthog::timer timer;


        // Pre-initialised variables to use in search().
        Successor* search_successors;
        SearchNode* search_nodes_to_push;

        void init()
        {
            verbose = false;
            search_successors = new Successor [mesh->max_poly_sides + 2];
            search_nodes_to_push = new SearchNode [mesh->max_poly_sides + 2];
            node_pool = new warthog::mem::cpool(sizeof(SearchNode));
            init_root_pruning();
        }
        void init_root_pruning()
        {
            assert(mesh != nullptr);
            search_id = 0;
            size_t num_vertices = mesh->mesh_vertices.size();
            root_g_values.resize(num_vertices);
            root_search_ids.resize(num_vertices);
            fill(root_search_ids.begin(), root_search_ids.end(), 0);
        }
        void init_search()
        {
            assert(node_pool);
            node_pool->reclaim();
            search_id++;
            open_list = pq();
            final_node = nullptr;
            nodes_generated = 0;
            nodes_pushed = 0;
            nodes_popped = 0;
            nodes_pruned_post_pop = 0;
            successor_calls = 0;
            insertion_time = 0;
            set_end_polygon();
            gen_initial_nodes();
        }

        void init_search(PointLocation start_pl ,PointLocation end_pl)
        {
            assert(node_pool);
            node_pool->reclaim();
            search_id++;
            open_list = pq();
            final_node = nullptr;
            nodes_generated = 0;
            nodes_pushed = 0;
            nodes_popped = 0;
            nodes_pruned_post_pop = 0;
            successor_calls = 0;

            set_end_polygon(end_pl);
            gen_initial_nodes(start_pl );
        }
        void set_end_polygon();
        void gen_initial_nodes();
        int succ_to_node(
                SearchNodePtr parent, Successor* successors,
                int num_succ, SearchNode* nodes
        );
        void print_node(SearchNodePtr node, std::ostream& outfile);

    public:
        int nodes_generated;        // Nodes stored in memory
        int nodes_pushed;           // Nodes pushed onto open
        int nodes_popped;           // Nodes popped off open
        int nodes_pruned_post_pop;  // Nodes we prune right after popping off
        int successor_calls;        // Times we call get_successors
        bool verbose;
        double sort_cost;

        double insertion_time;

        SearchInstance() = default;
        SearchInstance(MeshPtr m) : mesh(m) { init(); }
        SearchInstance(MeshPtr m, Point s, Point g) :
                mesh(m), start(s), goal(g) { init(); }
        SearchInstance(SearchInstance const &) = delete;
        void operator=(SearchInstance const &x) = delete;
        ~SearchInstance()
        {
            if (node_pool)
            {
                delete node_pool;
            }
            delete[] search_successors;
            delete[] search_nodes_to_push;
        }

        void set_start_goal(Point s, Point g)
        {
            start = s;
            goal = g;
            final_node = nullptr;
        }

        bool search();
        double get_cost()
        {
            if (final_node == nullptr)
            {
                return INF;
            }

            return final_node->f;
        }

        double get_search_micro()
        {
            return timer.elapsed_time_micro();
        }

        double get_search_nano()
        {
//            std::cout<<timer.elapsed_time_nano()<<std::endl;
//            std::cout<<timer.elapsed_time_micro()<<std::endl;
            return timer.elapsed_time_nano();
        }

        void get_path_points(std::vector<Point>& out);
        void print_search_nodes(std::ostream& outfile);


        bool search(PointLocation start_pl, PointLocation end_pl);

        void set_end_polygon(const PointLocation& end_pl);

        void gen_initial_nodes(const PointLocation& start_pl);

        bool bounded_search(PointLocation start_pl, PointLocation end_pl, double bound);

        void get_path_vertex_id(std::vector<int> &out);

        void get_path_vertex_id(std::vector<int> &out, int s_id, int g_id);

        bool search_covisible(PointLocation start_pl, PointLocation end_pl);

        bool triangle_search();

        bool triangle_search(vector<tuple<int, int>> &cur_triangle);
    };

}
