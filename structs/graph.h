#pragma once
#include<vector>
#include<fstream>
#include <algorithm>
#include<cstdint>
#include<iostream>
#include<assert.h>
#include<map>
#include <cmath>
#include <sstream>
#include <vec_io.h>
#include <random>
#include "geometry.h"

using namespace std;
namespace polyanya {

    // weighted directed Graph.
    class Graph {
    public:

        // sorted different value into vector for speed up.


        vector<int> vertices;
        vector<int> out_vertices;
        vector<double> distance_cost;


        int number_of_vertices;
        int number_of_edges;

        Graph() {
            number_of_vertices = 0;
            number_of_edges = 0;
            vertices.clear();
            out_vertices.clear();
            distance_cost.clear();
        }

        ~Graph() {
            vertices.clear();
            out_vertices.clear();
            distance_cost.clear();
        }

        void load_graph(const string &file_name);

        vector<int> generate_DFS_ordering();

        void resort_graph(const vector<int> &ordering);

    private:
        // for generate query only :
        bool is_connected(int s, int t);
    };
}
