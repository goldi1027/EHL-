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

namespace polyanya
{

// A search node.
// Only makes sense given a mesh and an endpoint, which the node does not store.
// This means that the f value needs to be set manually.
struct SearchNode
{
    SearchNode* parent;
    // Note that all Points here will be in terms of a Cartesian plane.
    int root; // -1 if start

    // If possible, set the orientation of left / root / right to be
    // "if I'm standing at 'root' and look at 'left', 'right' is on my right"
    Point left, right;

    // The left vertex of the edge the interval is lying on.
    // When generating the successors of this node, end there.
    int left_vertex;

    // The right vertex of the edge the interval is lying on.
    // When generating the successors of this node, start there.
    int right_vertex;

    // Index of the polygon we're going to "push" into.
    // Every successor must lie within this polygon.
    int next_polygon;

    double f, g;

    bool reached = false;
    int goal_id = -1;
    int heuristic_gid = -1;
    int edge_id = -1;



    SearchNode() {}
    SearchNode(SearchNode* p, int rid, Point l, Point r, int lv, int rv, int next_poly, double f, double g):
      parent(p), root(rid), left(l), right(r), left_vertex(lv), right_vertex(rv), next_polygon(next_poly), f(f), g(g) { }

    // Comparison.
    // Always take the "smallest" search node in a priority queue.
    bool operator<(const SearchNode& other) const
    {
        if (this->f == other.f)
        {
            // If two nodes have the same f, the one with the bigger g
            // is "smaller" to us.
            return this->g > other.g;
        }
        return this->f < other.f;
    }

    bool operator>(const SearchNode& other) const
    {
        if (this->f == other.f)
        {
            return this->g < other.g;
        }
        return this->f > other.f;
    }

    friend std::ostream& operator<<(std::ostream& stream, const SearchNode& sn)
    {
        return stream << "SearchNode ([" << sn.root << ", [" << sn.left << ", "
                      << sn.right << "]], f=" << sn.f << ", g=" << sn.g
                      << ", poly=" << sn.next_polygon << ")";
    }

    void set_reached() { reached = true; }
    void set_goal_id(int gid) { goal_id = gid; }
};

typedef SearchNode* SearchNodePtr;

}
