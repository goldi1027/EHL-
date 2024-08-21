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

using namespace std;
namespace polyanya
{

    template<typename T, typename Compare = std::greater<T> >
    struct PointerComp2
    {
        bool operator()(const T* x,
                        const T* y) const
        {
            return Compare()(*x, *y);
        }
    };

    typedef Mesh* MeshPtr;

// Polyanya instance for point to point search
    class visibleSearchInstance
    {
        typedef std::priority_queue<SearchNodePtr, std::vector<SearchNodePtr>,
                PointerComp2<SearchNode> > pq;
    private:
        warthog::mem::cpool* node_pool;
        MeshPtr mesh;




        // set by init_search
        pq open_list;
        int size;
        // Best g value for a specific vertex.
        std::vector<double> root_g_values;
        // Contains the current search id if the root has been reached by
        // the search.
        std::vector<int> root_search_ids;  // also used for root-level pruning

        std::vector<double> vertices_h_values;
        // Contains the current search id if the root has been reached by
        // the search.
        std::vector<int> vertices_search_ids;

        std::vector<int> visible_vertices_list;
        int vertice_search_id;

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
            vertice_search_id =0;

            size_t num_vertices = mesh->mesh_vertices.size();
            root_g_values.resize(num_vertices);
            root_search_ids.resize(num_vertices);
            fill(root_search_ids.begin(), root_search_ids.end(), 0);
            vertices_h_values.resize(num_vertices);
            vertices_search_ids.resize(num_vertices);
            fill(vertices_search_ids.begin(), vertices_search_ids.end(), 0);

            visible_vertices_list.resize(num_vertices);
            fill(visible_vertices_list.begin(), visible_vertices_list.end(), 0);

            size = num_vertices;
        }

        void init_search(PointLocation start_pl,vector<int>& visible_vertex)
        {
            assert(node_pool);
            node_pool->reclaim();
            search_id++;
            vertice_search_id ++;
            open_list = pq();
//            final_node = nullptr;
            nodes_generated = 0;
            nodes_pushed = 0;
            nodes_popped = 0;
            successor_calls = 0;
            gen_initial_nodes_for_visiblity_search(start_pl,visible_vertex);
        }


        int succ_to_node(
                SearchNodePtr parent, Successor* successors,
                int num_succ, SearchNode* nodes
        );
        void print_node(SearchNodePtr node, std::ostream& outfile);

    public:
        int degree;
        int nodes_generated;        // Nodes stored in memory
        int nodes_pushed;           // Nodes pushed onto open
        int nodes_popped;           // Nodes popped off open
        int successor_calls;        // Times we call get_successors
        bool verbose;
        int search_id;

        Point start, goal;

        visibleSearchInstance() = default;


        visibleSearchInstance(MeshPtr m ) : mesh(m) { init(); }

        visibleSearchInstance(MeshPtr m, Point s, Point g) :
                mesh(m), start(s), goal(g) { init(); }
        visibleSearchInstance(visibleSearchInstance const &) = delete;
        void operator=(visibleSearchInstance const &x) = delete;
        ~visibleSearchInstance()
        {
            if (node_pool)
            {
                delete node_pool;
            }
            delete[] search_successors;
            delete[] search_nodes_to_push;
        }

        void gen_initial_nodes_for_visiblity_search(const PointLocation &start_pl, vector<int> &visible_vertex);

        bool search_visible_vertices(int s, PointLocation start_pl, vector<int> &visible_vertex);


    };

}