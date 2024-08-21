/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#pragma once
#include <vector>
#include "edge.h"

namespace polyanya
{

struct Polygon
{
    // "int" here means an array index.
    std::vector<int> vertices;
    std::vector<int> polygons;
    std::vector<int> edges;
    // leave it later, we can use it to check the closest centorid

    bool is_one_way;
    double min_x, max_x, min_y, max_y;

    Edge getEdge(int index) const {
        Edge e;
        if(index == 0){
            e.vertices = std::make_pair(vertices[index],vertices[vertices.size()-1]);
        }else{
            e.vertices = std::make_pair (vertices[index],vertices[index-1]);
        }
        return  e;
    }
};

}
