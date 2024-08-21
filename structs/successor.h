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

struct Successor
{
    enum Type
    {
        RIGHT_NON_OBSERVABLE,
        OBSERVABLE,
        LEFT_NON_OBSERVABLE,
    };

    Type type;

    Point left, right;

    // Used to get next_polygon (Polygon.polygons[left_ind])
    // as well as right_vertex (Polygon.vertices[left_ind - 1])
    int poly_left_ind;

    friend std::ostream& operator<<(std::ostream& stream,
                                    const Successor& succ)
    {
        const std::string lookup[] = {
            "RIGHT_NON_OBSERVABLE",
            "OBSERVABLE",
            "LEFT_NON_OBSERVABLE",
        };
        return stream << "Successor(" << lookup[static_cast<int>(succ.type)]
                      << ", " << succ.left << " " << succ.right << ", "
                      << "poly_left=" << succ.poly_left_ind << ")";
    }
};

}
