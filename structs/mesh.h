/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#pragma once
#include "polygon.h"
#include "point.h"
#include "vertex.h"
#include "edge.h"
#include "successor.h"
#include "graph.h"
#include <vector>
#include <iostream>
#include <map>
#include <memory>



namespace polyanya
{

    struct Grid {
        std::map<int, std::vector<int>> centroid_map;
    };
    struct PolyContainment
    {
        enum Type
        {
            // Does not use any ints.
                    OUTSIDE,

            // Does not use any ints.
                    INSIDE,

            // Uses adjacent_poly, vertex1 and vertex2.
                    ON_EDGE,

            // Uses vertex1.
                    ON_VERTEX,
        };

        Type type;

        int adjacent_poly;

        // If on edge, vertex1/vertex2 represents the left/right vertices of the
        // edge when looking from a point in the poly.
        int vertex1, vertex2;

        friend std::ostream& operator<<(std::ostream& stream,
                                        const PolyContainment& pc)
        {
            switch (pc.type)
            {
                case PolyContainment::OUTSIDE:
                    return stream << "OUTSIDE";

                case PolyContainment::INSIDE:
                    return stream << "INSIDE";

                case PolyContainment::ON_EDGE:
                    return stream << "ON_EDGE (poly " << pc.adjacent_poly
                                  << ", vertices " << pc.vertex1 << ", "
                                  << pc.vertex2 << ")";

                case PolyContainment::ON_VERTEX:
                    return stream << "ON_VERTEX (" << pc.vertex1 << ")";

                default:
                    assert(false);
                    return stream;
            }
        }
    };


    struct PointLocation
    {
        enum Type
        {
            // Does not use any ints.
                    NOT_ON_MESH,

            // Uses poly1 (the polygon it is on).
                    IN_POLYGON,

            // Uses poly1 (the polygon it is on) and both vertices.
                    ON_MESH_BORDER,       // edge: a polygon is not traversable

            // Uses poly1, poly2 and both vertices.
                    ON_EDGE,              // edge: both polygons are traversable

            // Uses vertex1.
            // Can use poly1 to specify the "grid corrected poly".
            // Will need to manually assign poly1, though.
                    ON_CORNER_VERTEX_AMBIG,   // vertex; two+ polygons are not traversable

            // Uses vertex1. Also returns an arbirary traversable adjacent
            // polygon in poly1.
                    ON_CORNER_VERTEX_UNAMBIG, // vertex; one polygon is not traversable

            // Uses vertex1. Also returns an arbitrary adjacent polygon in poly1.
                    ON_NON_CORNER_VERTEX, // vertex: all polygons are traversable
        };

        Type type;
        int poly1, poly2;
        // If on edge, vertex1/vertex2 represents the left/right vertices of the
        // edge when looking from a point in poly1.
        int vertex1, vertex2;

        friend std::ostream& operator<<(std::ostream& stream,
                                        const PointLocation& pl)
        {
            switch (pl.type)
            {
                case PointLocation::NOT_ON_MESH:
                    return stream << "NOT_ON_MESH";

                case PointLocation::IN_POLYGON:
                    return stream << "IN_POLYGON (" << pl.poly1 << ")";

                case PointLocation::ON_MESH_BORDER:
                    return stream << "ON_MESH_BORDER (poly " << pl.poly1
                                  << ", vertices " << pl.vertex1 << ", "
                                  << pl.vertex2 << ")";

                case PointLocation::ON_EDGE:
                    return stream << "ON_EDGE (polys "
                                  << pl.poly1 << ", " << pl.poly2 << ", vertices "
                                  << pl.vertex1 << ", " << pl.vertex2 << ")";

                case PointLocation::ON_CORNER_VERTEX_AMBIG:
                    return stream << "ON_CORNER_VERTEX_AMBIG (" << pl.vertex1
                                  << ", poly? " << pl.poly1 << ")";

                case PointLocation::ON_CORNER_VERTEX_UNAMBIG:
                    return stream << "ON_CORNER_VERTEX_UNAMBIG (" << pl.vertex1
                                  << ", poly " << pl.poly1 << ")";

                case PointLocation::ON_NON_CORNER_VERTEX:
                    return stream << "ON_NON_CORNER_VERTEX (" << pl.vertex1
                                  << ", poly " << pl.poly1 << ")";

                default:
                    assert(false);
                    return stream;
            }
        }

        bool operator==(const PointLocation& other) const
        {
            if (type != other.type)
            {
                return false;
            }

            switch (type)
            {
                case PointLocation::NOT_ON_MESH:
                    return true;

                case PointLocation::IN_POLYGON:
                    return poly1 == other.poly1;

                case PointLocation::ON_MESH_BORDER:
                    if (poly1 != other.poly1)
                    {
                        return false;
                    }
                    if (vertex1 == other.vertex1 && vertex2 == other.vertex2)
                    {
                        return true;
                    }
                    if (vertex1 == other.vertex2 && vertex2 == other.vertex1)
                    {
                        return true;
                    }
                    return false;

                case PointLocation::ON_EDGE:
                    if (poly1 == other.poly1 && poly2 == other.poly2 &&
                        vertex1 == other.vertex1 && vertex2 == other.vertex2)
                    {
                        return true;
                    }
                    if (poly1 == other.poly2 && poly2 == other.poly1 &&
                        vertex1 == other.vertex2 && vertex2 == other.vertex1)
                    {
                        return true;
                    }
                    return false;

                case PointLocation::ON_CORNER_VERTEX_AMBIG:
                case PointLocation::ON_CORNER_VERTEX_UNAMBIG:
                case PointLocation::ON_NON_CORNER_VERTEX:
                    return vertex1 == other.vertex1;

                default:
                    assert(false);
                    return false;
            }
        }

        bool operator!=(const PointLocation& other) const
        {
            return !((*this) == other);
        }
    };

    class Mesh
    {
    private:
        std::map<double, std::vector<int>> slabs;
        double min_x, max_x, min_y, max_y;

    public:
        Mesh() { }
        Mesh(std::istream& infile);
        Mesh(std::istream& infile, int poly);
        std::vector<Vertex> mesh_vertices;
        std::vector<Polygon> mesh_polygons;
        std::vector<Edge> mesh_edges;
        std::vector<Edge> obstacle_edges;
        std::vector<int> cpd_to_vertices_mapper;
        std::vector<std::vector<int>> cpd_out_vertices;
        std::vector<std::vector<double>> cpd_distance_cost;

        int max_poly_sides;
        int width;
        int height;

        void read(std::istream& infile);
        void read_polygon(std::istream& infile);
        void precalc_point_location();
        void print(std::ostream& outfile);
        PolyContainment poly_contains_point(int poly, Point& p);
        PointLocation get_point_location(Point& p);
        PointLocation get_point_location_naive(Point& p);

        void print_polygon(std::ostream& outfile, int index);
        void print_vertex(std::ostream& outfile, int index);
        void print_edge(std::ostream& outfile, int index);
        void pre_calculate_edge();
        void Initialize_edge();
        void Assign_edge_ID(int polygon, Edge& e, int edgeID);
        double get_minx() { return min_x; }
        double get_maxx() { return max_x; }
        double get_miny() { return min_y; }
        double get_maxy() { return max_y; }


        void pre_compute_obstacle_edge_on_vertex();


        void mark_turning_point(const char *fname2);

        void mark_turning_point_polygon();
        void get_grid_width_height(int &map_width, int &map_height);

        void create_cpd_to_vertices_mapper(std::vector<int> vertices_mapper, std::vector<int> cpd_mapper);

        int check_edge_vertex(int vertex_1, int vertex_2, int polygon_id, int current_vertex);
        void create_cpd_out_list(polyanya::Graph g, const std::vector<int> &cpd_mapper);

        void get_polygon_width_height();

        int get_map_height();

        int get_map_width();

        double get_angles_between_points(const Point& p1, const Point& p2, const Point& p3);
    };

}
