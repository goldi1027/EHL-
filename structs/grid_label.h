/*
 * Implementation of EHL structure
 * Defines data structure of EHL
 */

#include <point.h>
#ifndef GRID_LABEL.H
#define GRID_LABEL.H



namespace polyanya{

    struct Compressed_convex_vertices_label{
        //compress convex_vertex, predecessor,
        //bool visilbility;
        unsigned compressed_label;
        double distance;
    };

    struct Convex_vertices_label{
        int convex_vertex;
        double distance;
        bool visibility;
        int predecessor;
    };
    struct Hub_label{
        int hub_id;
        vector<Convex_vertices_label> convex_labels;
        double min_lower_bound;
    };

    struct Grid_label{
        //four corner points
        Point a_p;
        Point b_p;
        Point c_p;
        Point d_p;
        vector<Hub_label> hub_labels;

    };


}

#endif GRID_LABEL.H