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
#include <iostream>
#include <vector>
#include <string>

namespace polyanya
{

struct Scenario
{
    int bucket;
    double xsize, ysize;
    Point start, goal;
    double gridcost;
};

bool load_scenarios(std::istream& infile, std::vector<Scenario>& out);

bool load_workload_scenarios(std::istream& infile, std::vector<Scenario>& out);

}
