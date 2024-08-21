/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#pragma once
#include "point.h"
#include "mesh.h"
#include <vector>

namespace polyanya
{

// A point in the polygon mesh.
struct Vertex
{
    Point p;
    // "int" here means an array index.
    std::vector<int> polygons;
    bool is_corner;
    bool is_ambig;
    std::vector<int> obstacle_edge;
    // since all the obstacle edges connect to the vertex itself, we only record the id of connected vertex;
    // obstacle edge : (this,Vertex[i]);
    bool is_turning_vertex = false;
    int cpd_id;

};

}
