/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#include "mesh.h"
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "geometry.h"
#include <fstream>
#include <sstream>
#include <grid_label.h>
#include "graph.h"

namespace polyanya
{


Mesh::Mesh(std::istream& infile)
{
    read(infile);
    precalc_point_location();
    Initialize_edge();
}

Mesh::Mesh(std::istream& infile, int poly)
{
    read_polygon(infile);
    precalc_point_location();
    Initialize_edge();
}

void Mesh::read(std::istream& infile)
{
    #define fail(message) std::cerr << message << std::endl; exit(1);
    std::string header;
    int version;

    if (!(infile >> header))
    {
        fail("Error reading header");
    }
    if (header != "mesh")
    {
        std::cerr << "Got header '" << header << "'" << std::endl;
        fail("Invalid header (expecting 'mesh')");
    }

    if (!(infile >> version))
    {
        fail("Error getting version number");
    }
    if (version != 2)
    {
        std::cerr << "Got file with version " << version << std::endl;
        fail("Invalid version (expecting 2)");
    }

    int V, P;
    if (!(infile >> V >> P))
    {
        fail("Error getting V and P");
    }
    if (V < 1)
    {
        std::cerr << "Got " << V << " vertices" << std::endl;
        fail("Invalid number of vertices");
    }
    if (P < 1)
    {
        std::cerr << "Got " << P << " polygons" << std::endl;
        fail("Invalid number of polygons");
    }

    mesh_vertices.resize(V);
    mesh_polygons.resize(P);


    for (int i = 0; i < V; i++)
    {
        Vertex& v = mesh_vertices[i];
        v.is_corner = false;
        v.is_ambig = false;
        if (!(infile >> v.p.x >> v.p.y))
        {
            fail("Error getting vertex point");
        }
        int neighbours;
        if (!(infile >> neighbours))
        {
            fail("Error getting vertex neighbours");
        }
        if (neighbours < 2)
        {
            std::cerr << "Got " << neighbours << " neighbours" << std::endl;
            fail("Invalid number of neighbours around a point");
        }
        v.polygons.resize(neighbours);
        for (int j = 0; j < neighbours; j++)
        {
            int polygon_index;
            if (!(infile >> polygon_index))
            {
                fail("Error getting a vertex's neighbouring polygon");
            }
            if (polygon_index >= P)
            {
                std::cerr << "Got a polygon index of " \
                          << polygon_index << std::endl;
                fail("Invalid polygon index when getting vertex");
            }
            v.polygons[j] = polygon_index;
            if (polygon_index == -1)
            {
                if (v.is_corner)
                {
                    if (!v.is_ambig)
                    {
                        v.is_ambig = true;
                    }
                }
                else
                {
                    v.is_corner = true;
                }
            }
        }
    }


    max_poly_sides = 0;
    for (int i = 0; i < P; i++)
    {
        Polygon& p = mesh_polygons[i];
        int n;
        if (!(infile >> n))
        {
            fail("Error getting number of vertices of polygon");
        }
        if (n < 3)
        {
            std::cerr << "Got " << n << " vertices" << std::endl;
            fail("Invalid number of vertices in polygon");
        }
        p.vertices.resize(n);
        p.polygons.resize(n);
        if (n > max_poly_sides)
        {
            max_poly_sides = n;
        }

        for (int j = 0; j < n; j++)
        {
            int vertex_index;
            if (!(infile >> vertex_index))
            {
                fail("Error getting a polygon's vertex");
            }
            if (vertex_index >= V)
            {
                std::cerr << "Got a vertex index of " \
                          << vertex_index << std::endl;
                fail("Invalid vertex index when getting polygon");
            }
            p.vertices[j] = vertex_index;
            if (j == 0)
            {
                p.min_x = mesh_vertices[vertex_index].p.x;
                p.min_y = mesh_vertices[vertex_index].p.y;
                p.max_x = mesh_vertices[vertex_index].p.x;
                p.max_y = mesh_vertices[vertex_index].p.y;
            }
            else
            {
                p.min_x = std::min(p.min_x, mesh_vertices[vertex_index].p.x);
                p.min_y = std::min(p.min_y, mesh_vertices[vertex_index].p.y);
                p.max_x = std::max(p.max_x, mesh_vertices[vertex_index].p.x);
                p.max_y = std::max(p.max_y, mesh_vertices[vertex_index].p.y);
            }
        }
        // mesh min/max
        if (i == 0)
        {
            min_x = p.min_x;
            min_y = p.min_y;
            max_x = p.max_x;
            max_y = p.max_y;
        }
        else
        {
            min_x = std::min(min_x, p.min_x);
            min_y = std::min(min_y, p.min_y);
            max_x = std::max(max_x, p.max_x);
            max_y = std::max(max_y, p.max_y);
        }

        bool found_trav = false;
        p.is_one_way = true;
        for (int j = 0; j < n; j++)
        {
            int polygon_index;
            if (!(infile >> polygon_index))
            {
                fail("Error getting a polygon's neighbouring polygon");
            }
            if (polygon_index >= P)
            {
                std::cerr << "Got a polygon index of " \
                          << polygon_index << std::endl;
                fail("Invalid polygon index when getting polygon");
            }
            if (polygon_index != -1)
            {
                if (found_trav)
                {
                    if (p.is_one_way)
                    {
                        p.is_one_way = false;
                    }
                }
                else
                {
                    found_trav = true;
                }
            }
            p.polygons[j] = polygon_index;
        }
    }

    double temp;
    if (infile >> temp)
    {
        fail("Error parsing mesh (read too much)");
    }
    #undef fail
}

void Mesh::read_polygon(std::istream& infile)
    {
#define fail(message) std::cerr << message << std::endl; exit(1);
        std::string header;
        int version;
        if (!(infile >> header))
        {
            fail("Error reading header");
        }
        if (header != "mesh")
        {
            std::cerr << "Got header '" << header << "'" << std::endl;
            fail("Invalid header (expecting 'mesh')");
        }

        if (!(infile >> version))
        {
            fail("Error getting version number");
        }
        if (version != 3)
        {
            std::cerr << "Got file with version " << version << std::endl;
            fail("Invalid version (expecting 2)");
        }

        int V, P;
        if (!(infile >> V >> P))
        {
            fail("Error getting V and P");
        }
        if (V < 1)
        {
            std::cerr << "Got " << V << " vertices" << std::endl;
            fail("Invalid number of vertices");
        }
        if (P < 1)
        {
            std::cerr << "Got " << P << " polygons" << std::endl;
            fail("Invalid number of polygons");
        }

        mesh_vertices.resize(V);
        mesh_polygons.resize(P);

        double cor_min_x = 0;
        double cor_min_y = 0;
        for (int i = 0; i < V; i++)
        {
            Vertex& v = mesh_vertices[i];
            v.is_corner = false;
            v.is_ambig = false;
            if (!(infile >> v.p.x >> v.p.y))
            {
                fail("Error getting vertex point");
            }
            if(v.p.x < cor_min_x && v.p.y < cor_min_y){
                cor_min_x = v.p.x;
                cor_min_y = v.p.y;
            }
        }

        max_poly_sides = 0;
        for (int i = 0; i < P; i++)
        {
            Polygon& p = mesh_polygons[i];
            int traversable;
            if(!(infile >> traversable)){
                fail("Error getting traversable info");
            }
            if(traversable == 0){
                int read_n;
                infile >> read_n;
                for(int j = 0; j < read_n; ++j){
                    int read_vertex;
                    infile >> read_vertex;
                    read_vertex = read_vertex - 1;
//                    mesh_vertices[read_vertex].polygons.push_back(-1);
                }
                for(int j = 0; j < read_n; ++j){
                    int read_poly;
                    infile >> read_poly;
                }
                continue;
            }
            int n;
            if (!(infile >> n))
            {
                fail("Error getting number of vertices of polygon");
            }
            if (n < 3)
            {
                std::cerr << "Got " << n << " vertices" << std::endl;
                //cout << i << endl;
                fail("Invalid number of vertices in polygon");
            }
            p.vertices.resize(n);
            p.polygons.resize(n);
//            if (n > max_poly_sides)
//            {
//                max_poly_sides = n;
//            }

            for (int j = 0; j < n; j++)
            {
                int vertex_index;
                if (!(infile >> vertex_index))
                {
                    fail("Error getting a polygon's vertex");
                }
                vertex_index = vertex_index - 1;
                if (vertex_index >= V)
                {
                    std::cerr << "Got a vertex index of " \
                          << vertex_index << std::endl;
                    fail("Invalid vertex index when getting polygon");
                }
                p.vertices[j] = vertex_index;
                if (j == 0)
                {
                    p.min_x = mesh_vertices[vertex_index].p.x;
                    p.min_y = mesh_vertices[vertex_index].p.y;
                    p.max_x = mesh_vertices[vertex_index].p.x;
                    p.max_y = mesh_vertices[vertex_index].p.y;
                }
                else
                {
                    p.min_x = std::min(p.min_x, mesh_vertices[vertex_index].p.x);
                    p.min_y = std::min(p.min_y, mesh_vertices[vertex_index].p.y);
                    p.max_x = std::max(p.max_x, mesh_vertices[vertex_index].p.x);
                    p.max_y = std::max(p.max_y, mesh_vertices[vertex_index].p.y);
                }
            }
            // mesh min/max
            if (i == 0)
            {
                min_x = p.min_x;
                min_y = p.min_y;
                max_x = p.max_x;
                max_y = p.max_y;
            }
            else
            {
                min_x = std::min(min_x, p.min_x);
                min_y = std::min(min_y, p.min_y);
                max_x = std::max(max_x, p.max_x);
                max_y = std::max(max_y, p.max_y);
            }

            bool found_trav = false;
            p.is_one_way = true;
            for (int j = 0; j < n; j++)
            {
                int polygon_index;
                if (!(infile >> polygon_index))
                {
                    fail("Error getting a polygon's neighbouring polygon");
                }
                if(polygon_index < 0){
                    polygon_index = -1;
                }
                else{
                    polygon_index = polygon_index - 1;
                }
                if (polygon_index >= P)
                {
                    std::cerr << "Got a polygon index of " \
                          << polygon_index << std::endl;
                    fail("Invalid polygon index when getting polygon");
                }
                if (polygon_index != -1)
                {
                    if (found_trav)
                    {
                        if (p.is_one_way)
                        {
                            p.is_one_way = false;
                        }
                    }
                    else
                    {
                        found_trav = true;
                    }
                }
                p.polygons[j] = polygon_index;
            }
        }

        double temp;
        if (infile >> temp)
        {
            fail("Error parsing mesh (read too much)");
        }

        for(int i = 0; i < mesh_vertices.size(); ++i){
            mesh_vertices[i].p.x = mesh_vertices[i].p.x + abs(cor_min_x);
            mesh_vertices[i].p.y = mesh_vertices[i].p.y + abs(cor_min_y);
        }
        for(int i = 0; i < mesh_polygons.size(); ++i){
            for(int j = 0; j < mesh_polygons[i].vertices.size(); ++j){
                int vertex_id = mesh_polygons[i].vertices[j];
                mesh_vertices[vertex_id].polygons.push_back(i);
            }
        }

        for(int i = 0; i < mesh_vertices.size(); ++i){
            map<int, int> neighbour_info;
            vector<pair<double,int>> temp_vertex;
            Point start_point;
            start_point.x = mesh_vertices[i].p.x + 1;
            start_point.y = mesh_vertices[i].p.y + 1;
            for(int j = 0; j < mesh_vertices[i].polygons.size(); ++j){
                int current_poly = mesh_vertices[i].polygons[j];
                for(int k = 0; k < mesh_polygons[current_poly].vertices.size(); ++k){
                    if(mesh_polygons[current_poly].vertices[k] == i){
                        if(k == 0){
//                            if(j == 0) {
//                                start_point = mesh_vertices[mesh_polygons[current_poly].vertices[
//                                        mesh_polygons[current_poly].vertices.size() - 1]].p;
//                            }
                            double angle = get_angles_between_points(mesh_vertices[i].p, start_point, mesh_vertices[mesh_polygons[current_poly].vertices[mesh_polygons[current_poly].vertices.size() -1]].p);
//                            temp_vertex.push_back(make_pair(angle, mesh_polygons[current_poly].vertices[
//                                        mesh_polygons[current_poly].vertices.size() - 1]));
                            temp_vertex.push_back(make_pair(angle, current_poly));
                            neighbour_info.insert({current_poly,mesh_polygons[current_poly].polygons[k]});
                        }
                        else{
                            double angle = get_angles_between_points(mesh_vertices[i].p, start_point,mesh_vertices[mesh_polygons[current_poly].vertices[k-1]].p);
                            temp_vertex.push_back(make_pair(angle,current_poly));
                            neighbour_info.insert({current_poly,mesh_polygons[current_poly].polygons[k]});
                        }
                    }
                }

            }
            vector<int> final_poly;
            sort(temp_vertex.begin(), temp_vertex.end());
            for(int j = 0; j < neighbour_info.size(); ++j){
                if(neighbour_info[temp_vertex[j].second] <= 0){
                    final_poly.push_back(-1);
                    final_poly.push_back(temp_vertex[j].second);
                }
                else{
                    final_poly.push_back(temp_vertex[j].second);
                }
            }
            mesh_vertices[i].polygons = final_poly;

        }

#undef fail
    }
    double Mesh::get_angles_between_points(const Point& p1, const Point& p2, const Point& p3){
        double AB = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
        double BC = sqrt(pow(p2.x - p3.x, 2) + pow(p2.y - p3.y, 2));
        double AC = sqrt(pow(p3.x - p1.x, 2) + pow(p3.y - p1.y, 2));
        double raidus = acos((BC * BC + AB * AB - AC * AC) / (2 * BC * AB));
        double degree = raidus * 180 / PI;
        if (get_orientation(p1, p2, p3) == Orientation::CCW) {
            degree = 360 - degree;
        }
        return 360 - degree;
    }

void Mesh::precalc_point_location()
{
    for (Vertex& v : mesh_vertices)
    {
        slabs[v.p.x] = std::vector<int>(0); // initialises the vector
    }
    for (int i = 0; i < (int) mesh_polygons.size(); i++)
    {
        const Polygon& p = mesh_polygons[i];
        const auto low_it = slabs.lower_bound(p.min_x);
        const auto high_it = slabs.upper_bound(p.max_x);

        for (auto it = low_it; it != high_it; it++)
        {
            it->second.push_back(i);
        }
    }
    for (auto& pair : slabs)
    {
        std::sort(pair.second.begin(), pair.second.end(),
            [&](const int& a, const int& b) -> bool
            {
                // Sorts based on the midpoints.
                // If tied, sort based on width of poly.
                const Polygon& ap = mesh_polygons[a], bp = mesh_polygons[b];
                const double as = ap.min_y + ap.max_y, bs = bp.min_y + bp.max_y;
                if (as == bs) {
                    return (ap.max_y - ap.min_y) > (bp.max_y - bp.min_y);
                }
                return as < bs;
            }
        );
    }
}

// Finds out whether the polygon specified by "poly" contains point P.
PolyContainment Mesh::poly_contains_point(int poly, Point& p)
{
    // The below is taken from
    // "An Efficient Test for a Point to Be in a Convex Polygon"
    // from the Wolfram Demonstrations Project
    // demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/

    // Assume points are in counterclockwise order.
    const Polygon& poly_ref = mesh_polygons[poly];
    if (p.x < poly_ref.min_x - EPSILON || p.x > poly_ref.max_x + EPSILON ||
        p.y < poly_ref.min_y - EPSILON || p.y > poly_ref.max_y + EPSILON)
    {
        return {PolyContainment::OUTSIDE, -1, -1, -1};
    }
    const Point& last_point_in_poly = mesh_vertices[poly_ref.vertices.back()].p;
    const Point ZERO = {0, 0};

    Point last = last_point_in_poly - p;
    if (last == ZERO)
    {
        return {PolyContainment::ON_VERTEX, -1, poly_ref.vertices.back(), -1};
    }

    int last_index = poly_ref.vertices.back();
    for (int i = 0; i < (int) poly_ref.vertices.size(); i++)
    {
        const int point_index = poly_ref.vertices[i];
        const Point cur = mesh_vertices[point_index].p - p;
        if (cur == ZERO)
        {
            return {PolyContainment::ON_VERTEX, -1, point_index, -1};
        }
        const double cur_a = last * cur;
        if (std::abs(cur_a) < EPSILON)
        {
            // The line going from cur to last goes through p.
            // This means that they are collinear.
            // The associated polygon should simply be polygons[i] in version
            // 2 of the file format.

            // Ensure that cur = c*last where c is negative.
            // If not, this means that the point is either outside or that this
            // segment is collinear to an adjacent one.
            if (cur.x)
            {
                if (!((cur.x > 0) ^ (last.x > 0)))
                {
                    last = cur;
                    last_index = point_index;
                    continue;
                }
            }
            else
            {
                if (!((cur.y > 0) ^ (last.y > 0)))
                {
                    last = cur;
                    last_index = point_index;
                    continue;
                }
            }
            return {PolyContainment::ON_EDGE, poly_ref.polygons[i],
                    point_index, last_index};
        }

        // Because we assume that the points are counterclockwise,
        // we can immediately terminate when we see a negatively signed area.
        if (cur_a < 0)
        {
            return {PolyContainment::OUTSIDE, -1, -1, -1};
        }
        last = cur;
        last_index = point_index;
    }
    return {PolyContainment::INSIDE, -1, -1, -1};
}

// Finds where the point P lies in the mesh.
PointLocation Mesh::get_point_location(Point& p)
{
    if (p.x < min_x - EPSILON || p.x > max_x + EPSILON ||
        p.y < min_y - EPSILON || p.y > max_y + EPSILON)
    {
        return {PointLocation::NOT_ON_MESH, -1, -1, -1, -1};
    }
    auto slab = slabs.upper_bound(p.x);
    if (slab == slabs.begin())
    {
        return {PointLocation::NOT_ON_MESH, -1, -1, -1, -1};
    }
    slab--;
    const std::vector<int>& polys = slab->second;
    const auto close_it = std::lower_bound(polys.begin(), polys.end(), p.y,
        [&](const int& poly_index, const double& y_coord) -> bool
        {
            // Sorts based on the midpoints.
            // If tied, sort based on width of poly.
            const Polygon& poly = mesh_polygons[poly_index];
            return poly.min_y + poly.max_y < y_coord * 2;
        }
    );
    const int close_index = close_it - polys.begin()
                            - (close_it == polys.end());
    // The plan is to take an index and repeatedly do:
    // +1, -2, +3, -4, +5, -6, +7, -8, ...
    // until it hits the edge. If it hits an edge, instead iterate normally.
    const int ps = polys.size();
    int i = close_index;
    int next_delta = 1;
    int walk_delta = 0; // way to go when walking normally

    while (i >= 0 && i < ps)
    {
        const int polygon = polys[i];
        const PolyContainment result = poly_contains_point(polygon, p);
        switch (result.type)
        {
            case PolyContainment::OUTSIDE:
                // Does not contain: try the next one.
                break;

            case PolyContainment::INSIDE:
                // This one strictly contains the point.
                return {PointLocation::IN_POLYGON, polygon, -1, -1, -1};

            case PolyContainment::ON_EDGE:
                // This one lies on the edge.
                // Chek whether the other one is -1.
                return {
                    (result.adjacent_poly == -1 ?
                     PointLocation::ON_MESH_BORDER :
                     PointLocation::ON_EDGE),
                    polygon, result.adjacent_poly,
                    result.vertex1, result.vertex2
                };

            case PolyContainment::ON_VERTEX:
                // This one lies on a corner.
            {

                const Vertex& v = mesh_vertices[result.vertex1];
                if (v.is_corner)
                {
                    if (v.is_ambig)
                    {
                        return {PointLocation::ON_CORNER_VERTEX_AMBIG, -1, -1,
                                result.vertex1, -1};
                    }
                    else
                    {
                        return {PointLocation::ON_CORNER_VERTEX_UNAMBIG,
                                polygon, -1, result.vertex1, -1};
                    }
                }
                else
                {
                    return {PointLocation::ON_NON_CORNER_VERTEX,
                            polygon, -1,
                            result.vertex1, -1};
                }
            }

            default:
                // This should not be reachable
                assert(false);
        }


        // do stuff
        if (walk_delta == 0)
        {
            const int next_i = i + next_delta * (2 * (next_delta & 1) - 1);
            if (next_i < 0)
            {
                // was going to go too far to the left.
                // start going right
                walk_delta = 1;
            }
            else if (next_i >= ps)
            {
                walk_delta = -1;
            }
            else
            {
                i = next_i;
                next_delta++;
            }
        }

        if (walk_delta != 0)
        {
            i += walk_delta;
        }
    }
    // Haven't returned yet, therefore P does not lie on the mesh.
    return {PointLocation::NOT_ON_MESH, -1, -1, -1, -1};
}

PointLocation Mesh::get_point_location_naive(Point& p)
{
    for (int polygon = 0; polygon < (int) mesh_polygons.size(); polygon++)
    {
        const PolyContainment result = poly_contains_point(polygon, p);
        switch (result.type)
        {
            case PolyContainment::OUTSIDE:
                // Does not contain: try the next one.
                break;

            case PolyContainment::INSIDE:
                // This one strictly contains the point.
                return {PointLocation::IN_POLYGON, polygon, -1, -1, -1};

            case PolyContainment::ON_EDGE:
                // This one lies on the edge.
                // Chek whether the other one is -1.
                return {
                    (result.adjacent_poly == -1 ?
                     PointLocation::ON_MESH_BORDER :
                     PointLocation::ON_EDGE),
                    polygon, result.adjacent_poly,
                    result.vertex1, result.vertex2
                };

            case PolyContainment::ON_VERTEX:
                // This one lies on a corner.
            {
                const Vertex& v = mesh_vertices[result.vertex1];
                if (v.is_corner)
                {
                    if (v.is_ambig)
                    {
                        return {PointLocation::ON_CORNER_VERTEX_AMBIG, -1, -1,
                                result.vertex1, -1};
                    }
                    else
                    {
                        return {PointLocation::ON_CORNER_VERTEX_UNAMBIG,
                                polygon, -1, result.vertex1, -1};
                    }
                }
                else
                {
                    return {PointLocation::ON_NON_CORNER_VERTEX,
                            polygon, -1,
                            result.vertex1, -1};
                }
            }

            default:
                // This should not be reachable
                assert(false);
        }
    }
    // Haven't returned yet, therefore P does not lie on the mesh.
    return {PointLocation::NOT_ON_MESH, -1, -1, -1, -1};
}

void Mesh::print(std::ostream& outfile)
{
    outfile << "mesh with " << mesh_vertices.size() << " vertices, " \
            << mesh_polygons.size() << " polygons" << std::endl;
    outfile << "vertices:" << std::endl;
    for (Vertex vertex : mesh_vertices)
    {
        outfile << vertex.p << " " << vertex.is_corner <<" " << vertex.is_ambig << std::endl;
        for(int i : vertex.polygons ){
            print_polygon(outfile,i);
            outfile<< i <<" " << std::endl;
        }

    }
    outfile << std::endl;
    outfile << "polygons:" << std::endl;
    for (Polygon polygon : mesh_polygons)
    {
        for (int vertex : polygon.vertices)
        {
            outfile << mesh_vertices[vertex].p << " ";
        }
        for (int polygon : polygon.polygons)
        {
            outfile << polygon << " " ;
        }
        outfile << std::endl;
        int index = 0 ;
        for (int edge : polygon.edges){
            if(edge == -1 ){
                Edge e = polygon.getEdge(index);
                outfile << mesh_vertices[e.vertices.first].p << mesh_vertices[e.vertices.second].p << std::endl;
            }else {
                Edge e = mesh_edges[edge];
                outfile << mesh_vertices[e.vertices.first].p << mesh_vertices[e.vertices.second].p << std::endl;
            }
            index ++;
        }
    }
    outfile << "Obstacle Edges"<< std::endl;
    for( auto edges : obstacle_edges){
        outfile << mesh_vertices[edges.vertices.first].p << mesh_vertices[edges.vertices.second].p <<std::endl;
    }

    outfile << "Mesh Edges"<< std::endl;
    for( auto edges : mesh_edges){
        outfile << mesh_vertices[edges.vertices.first].p << mesh_vertices[edges.vertices.second].p <<std::endl;
    }
}

void Mesh::print_polygon(std::ostream& outfile, int index)
{
    if (index == -1)
    {
        outfile << "P!!!";
        return;
    }
    outfile << "P" << index << " [";
    Polygon& poly = mesh_polygons[index];
    const auto vertices = poly.vertices;
    const int size = (int) vertices.size();
    for (int i = 0; i < size; i++)
    {
        print_vertex(outfile, vertices[i]);
        if (i != size - 1)
        {
            outfile << ", ";
        }
    }
    outfile << "]";
}

void Mesh::print_vertex(std::ostream& outfile, int index)
{
    outfile << "V" << index << " " <<  mesh_vertices[index].p<<std::endl;
    for(int o :  mesh_vertices[index].obstacle_edge){
        outfile << "obstacle_edg:" << mesh_vertices[index].p << " " <<  mesh_vertices[o].p<<std::endl;
    }
}

void Mesh::print_edge(std::ostream& outfile, int index)
{
    for( int i : mesh_polygons[index].edges){

        outfile << "E" << " " <<  i;
    }

    outfile << std::endl;
}



void Mesh::Initialize_edge(){
//    std::cout<< "initializing" << std::endl;
//    std::cout<< "Map Size"<< min_x << " "<< min_y << " "<< max_x <<" " << max_y << std::endl;
    for(int i =0; i < mesh_polygons.size(); i++){
        Polygon& p = mesh_polygons[i];
        p.edges.resize(p.vertices.size());
        //index record the vertex id,
        int obstacleEdgeId = 0;
        int index = 0;
        for(int polygon : p.polygons){
            if(polygon == -1){
                //obstacle here ;
                Edge obstacleEdges;
                if(index == 0){
                    obstacleEdges.vertices = std::make_pair (p.vertices[index],p.vertices[p.polygons.size()-1]);
                }else{
                    obstacleEdges.vertices = std::make_pair (p.vertices[index],p.vertices[index-1]);
                }
//                obstacleEdges.isObstacleEdge = true;
                p.edges[index] =  - 1;
                obstacle_edges.push_back(obstacleEdges);
                obstacleEdgeId ++;
            }else{
                // traversable edge but we leave it as  -1 * INF for now;
                p.edges[index] = -2;
            }
            index ++;
        }
    }
    pre_calculate_edge();
}


void Mesh::pre_calculate_edge(){
    // calculate traversable edge here;
    int edgeID = 0 ;
    // edge ID in mesh_edges;
    int polygonID = 0;
    for(Polygon& p : mesh_polygons){
        int index = 0;
        for(int edge: p.edges){
            if(edge ==  -2){
                //an edge that have not been assigned;
                Edge traverableEdges;
                if(index == 0){
                    traverableEdges.vertices = std::make_pair (p.vertices[index],p.vertices[p.vertices.size()-1]);
                }else{
                    traverableEdges.vertices = std::make_pair (p.vertices[index],p.vertices[index-1]);
                }
                p.edges[index] = edgeID;
                traverableEdges.polygons.first = polygonID;
                traverableEdges.polygons.second =p.polygons[index];
                //now assign the edge to adjacent polygon;
                Assign_edge_ID(p.polygons[index],traverableEdges,edgeID);
                mesh_edges.push_back(traverableEdges);
                edgeID ++;
            }
            index++;
        }
        polygonID++;

    }
}

void Mesh::Assign_edge_ID(int polygon, Edge& e, int edgeID){
    Polygon& checkPolygon = mesh_polygons[polygon];
    int size = checkPolygon.vertices.size();
    for(int i =0; i < size; i++){
        if(i==0){
            if (e.vertices.first == checkPolygon.vertices[0]  &&  e.vertices.second == checkPolygon.vertices[size-1]){
                checkPolygon.edges[0] = edgeID;
                break;
            }
            if (e.vertices.second == checkPolygon.vertices[0]  &&  e.vertices.first == checkPolygon.vertices[size-1]){
                checkPolygon.edges[0] = edgeID;
                break;
            }
        }else{
            if(e.vertices.first == checkPolygon.vertices[i] &&  e.vertices.second == checkPolygon.vertices[i-1]){
                checkPolygon.edges[i] = edgeID;
                break;
            }
            if(e.vertices.second == checkPolygon.vertices[i] &&  e.vertices.first == checkPolygon.vertices[i-1]){
                checkPolygon.edges[i] = edgeID;
                break;
            }

        }
    }
}



//    void Mesh::zordering(){
//        int i = 0;
//        for(Point& c : centroid_list){
//            uint_fast32_t x = c.x;
//            uint_fast32_t y = c.y;
//            uint_fast64_t zvalue = libmorton::morton2D_64_encode(x,y);
//            c.zvalue = zvalue;
//            c.centroidID = i;
//            //for debug delete it later;
//            zorder_list.push_back(c);
//            i++;
//        }
//        sort(zorder_list.begin( ), zorder_list.end( ), [ ]( const auto& lhs, const auto& rhs )
//        {
//            return lhs.zvalue < rhs.zvalue;
//        });
//        int id = 0;
//        for(Point& c : zorder_list){
//            centroid_list[c.centroidID].zorderID = id;
//            id++;
//        }
//    }

    void Mesh::pre_compute_obstacle_edge_on_vertex(){
        for(int i =0; i < mesh_polygons.size(); i++){
            Polygon& p = mesh_polygons[i];
//            p.edges.resize(p.vertices.size());
//            //index record the vertex id,
//            int obstacleEdgeId = 0;
            int index = 0;
            for(int polygon : p.polygons){
                if(polygon == -1){
                    //obstacle here ;
                    int v_id1;
                    int v_id2;
                    if(index == 0){
                        v_id1 = p.vertices[index];
                        v_id2 = p.vertices[p.polygons.size()-1];
                    }else{
                        v_id1 = p.vertices[index];
                        v_id2 = p.vertices[index-1];
                    }
//                obstacleEdges.isObstacleEdge = true;
                    mesh_vertices[v_id1].obstacle_edge.push_back(v_id2);
                    mesh_vertices[v_id2].obstacle_edge.push_back(v_id1);
                }
                index ++;
            }
        }
    }

    void Mesh::get_grid_width_height(int &map_width, int &map_height){
        map_width = width;
        map_height = height;

    }

    void Mesh::mark_turning_point_polygon(){
        int id = 0;
        for(Vertex& v: mesh_vertices){
            if(v.obstacle_edge.size()!=2 || v.is_ambig){
                id++;
                continue;
            }
            for(int i = 0; i < v.polygons.size(); ++i){
                int polygon_id = v.polygons[i];
                if(polygon_id == -1){
                    continue;
                }
                for(int j = 0; j < mesh_polygons[polygon_id].vertices.size(); ++j){
                    if(mesh_polygons[polygon_id].vertices[j] == id and mesh_polygons[polygon_id].edges[j] == -1){
                        if(j == mesh_polygons[polygon_id].vertices.size()-1){
                            if(mesh_polygons[polygon_id].vertices[j-1]== v.obstacle_edge[0] or mesh_polygons[polygon_id].vertices[j-1] == v.obstacle_edge[1]){
                                int obstacle_vertex;
                                if(mesh_polygons[polygon_id].vertices[j-1] == v.obstacle_edge[0]){
                                    obstacle_vertex = v.obstacle_edge[1];
                                }
                                else{
                                    obstacle_vertex = v.obstacle_edge[0];
                                }
                                int start_vertex = mesh_polygons[polygon_id].vertices[j-1];
                                int poly_vertex = mesh_polygons[polygon_id].vertices[0];
                                const Orientation& o1 = get_orientation(mesh_vertices[start_vertex].p, v.p, mesh_vertices[obstacle_vertex].p);
                                const Orientation& o2 = get_orientation(mesh_vertices[start_vertex].p,v.p, mesh_vertices[poly_vertex].p);
                                if(o1 != o2){
                                    v.is_turning_vertex = true;
                                }
                            }
                            else if(mesh_polygons[polygon_id].vertices[0]== v.obstacle_edge[0] or mesh_polygons[polygon_id].vertices[0] == v.obstacle_edge[1]){
                                int obstacle_vertex;
                                if(mesh_polygons[polygon_id].vertices[0] == v.obstacle_edge[0]){
                                    obstacle_vertex = v.obstacle_edge[1];
                                }
                                else{
                                    obstacle_vertex = v.obstacle_edge[0];
                                }
                                int start_vertex = mesh_polygons[polygon_id].vertices[0];
                                int poly_vertex = mesh_polygons[polygon_id].vertices[1];
                                const Orientation& o1 = get_orientation(mesh_vertices[start_vertex].p,v.p, mesh_vertices[obstacle_vertex].p);
                                const Orientation& o2 = get_orientation(mesh_vertices[start_vertex].p, v.p, mesh_vertices[poly_vertex].p);
                                if(o1 != o2){
                                    v.is_turning_vertex = true;
                                }
                            }
                        }
                        else if(j == 0){
                            if(mesh_polygons[polygon_id].vertices[j+1]== v.obstacle_edge[0] or mesh_polygons[polygon_id].vertices[j+1] == v.obstacle_edge[1]){
                                int obstacle_vertex;
                                if(mesh_polygons[polygon_id].vertices[j+1] == v.obstacle_edge[0]){
                                    obstacle_vertex = v.obstacle_edge[1];
                                }
                                else{
                                    obstacle_vertex = v.obstacle_edge[0];
                                }
                                int start_vertex = mesh_polygons[polygon_id].vertices[j+1];
                                int poly_vertex = mesh_polygons[polygon_id].vertices[2];
                                const Orientation& o1 = get_orientation(mesh_vertices[start_vertex].p, v.p, mesh_vertices[obstacle_vertex].p);
                                const Orientation& o2 = get_orientation(mesh_vertices[start_vertex].p, v.p,mesh_vertices[poly_vertex].p);
                                if(o1 != o2){
                                    v.is_turning_vertex = true;
                                }
                            }
                            else if(mesh_polygons[polygon_id].vertices[mesh_polygons[polygon_id].vertices.size()-1]== v.obstacle_edge[0] or mesh_polygons[polygon_id].vertices[mesh_polygons[polygon_id].vertices.size()-1] == v.obstacle_edge[1]){
                                int obstacle_vertex;
                                if(mesh_polygons[polygon_id].vertices[mesh_polygons[polygon_id].vertices.size()-1] == v.obstacle_edge[0]){
                                    obstacle_vertex = v.obstacle_edge[1];
                                }
                                else{
                                    obstacle_vertex = v.obstacle_edge[0];
                                }
                                int start_vertex = mesh_polygons[polygon_id].vertices[mesh_polygons[polygon_id].vertices.size()-1];
                                int poly_vertex = mesh_polygons[polygon_id].vertices[j+1];
                                const Orientation& o1 = get_orientation(mesh_vertices[start_vertex].p,v.p, mesh_vertices[obstacle_vertex].p);
                                const Orientation& o2 = get_orientation(mesh_vertices[start_vertex].p, v.p, mesh_vertices[poly_vertex].p);
                                if(o1 != o2){
                                    v.is_turning_vertex = true;
                                }
                            }
                        }
                        else{
                            if(mesh_polygons[polygon_id].vertices[j-1]== v.obstacle_edge[0] or mesh_polygons[polygon_id].vertices[j-1] == v.obstacle_edge[1]){
                                int obstacle_vertex;
                                if(mesh_polygons[polygon_id].vertices[j-1] == v.obstacle_edge[0]){
                                    obstacle_vertex = v.obstacle_edge[1];
                                }
                                else{
                                    obstacle_vertex = v.obstacle_edge[0];
                                }
                                int start_vertex = mesh_polygons[polygon_id].vertices[j-1];
                                int poly_vertex = mesh_polygons[polygon_id].vertices[j+1];
                                const Orientation& o1 = get_orientation(mesh_vertices[start_vertex].p, v.p, mesh_vertices[obstacle_vertex].p);
                                const Orientation& o2 = get_orientation(mesh_vertices[start_vertex].p,v.p,mesh_vertices[poly_vertex].p);
                                if(o1 != o2){
                                    v.is_turning_vertex = true;
                                }
                            }
                            else if(mesh_polygons[polygon_id].vertices[j+1]== v.obstacle_edge[0] or mesh_polygons[polygon_id].vertices[j+1] == v.obstacle_edge[j+1]){
                                int obstacle_vertex;
                                if(mesh_polygons[polygon_id].vertices[0] == v.obstacle_edge[0]){
                                    obstacle_vertex = v.obstacle_edge[1];
                                }
                                else{
                                    obstacle_vertex = v.obstacle_edge[0];
                                }
                                int start_vertex = mesh_polygons[polygon_id].vertices[j+1];
                                int poly_vertex = mesh_polygons[polygon_id].vertices[j-1];
                                const Orientation& o1 = get_orientation(mesh_vertices[start_vertex].p,v.p, mesh_vertices[obstacle_vertex].p);
                                const Orientation& o2 = get_orientation(mesh_vertices[start_vertex].p, v.p,mesh_vertices[poly_vertex].p);
                                if(o1 != o2){
                                    v.is_turning_vertex = true;
                                }
                            }
                        }
//                        if(o1 == o2){
//                            v.is_turning_vertex = true;
//                        }
                    }
                }
            }
            id++;
        }

}

    void Mesh::get_polygon_width_height() {
        int max_width = 0;
        int max_height = 0;
        for(Vertex& v: mesh_vertices){
            if(v.p.x > max_width){
                max_width = ceil(v.p.x);
            }
            if(v.p.y > max_height){
                max_height = ceil(v.p.y);
            }
        }
        width = max_width;
        height = max_height;
}

    void Mesh::mark_turning_point(const char* fname2){

        std::vector<std::pair<bool,int>> map;
        polyanya::LoadMap(fname2,map,width,height);

        int id = 0;
        for(Vertex& v :mesh_vertices){
            if(v.obstacle_edge.size()!=2 ||v.is_ambig){
                id ++;
                continue;
            }
            const Vertex& v1  = mesh_vertices[v.obstacle_edge[0]];
            const Vertex& v2  = mesh_vertices[v.obstacle_edge[1]];
            const Point middle = {(v1.p.x+v2.p.x)/2 , (v1.p.y+v2.p.y)/2  };
            double x = v.p.x;
            double y = v.p.y;

            if(middle.x> x){
                if(middle.y>y){
                    //up right
                    double key = y*width +x ;
                    if(!map[key].first){
                        mesh_vertices[id].is_turning_vertex = true;
                    }
                }else{
                    //bottom right
                    double key = (y-1)*width +x ;
                    if(!map[key].first){
                        mesh_vertices[id].is_turning_vertex = true;
                    }
                }

            }else{
                if(middle.y>y){
                    //up_left
                    double key = y*width +x-1 ;
                    if(!map[key].first){
                        mesh_vertices[id].is_turning_vertex = true;
                    }
                }else{
                    //bottom left
                    double key = (y-1) *width + x-1 ;
                    if(!map[key].first){
                        mesh_vertices[id].is_turning_vertex = true;
                    }
                }
            }
            id++;
        }
}    void Mesh::create_cpd_to_vertices_mapper(std::vector<int> vertices_mapper,std::vector<int> cpd_mapper) {
    // vertice mapper map vertices_id to visiblity graph_id
    // cpd mapper map graph_id to cpd_id

    cpd_to_vertices_mapper.clear();
    cpd_to_vertices_mapper.resize(cpd_mapper.size());
    for(int i = 0; i<vertices_mapper.size(); i++){
        if(vertices_mapper[i] == -1){
            mesh_vertices[i].cpd_id = -1;

        }else{
            // map vertice_id to cpd_id
            mesh_vertices[i].cpd_id = cpd_mapper[vertices_mapper[i]];
            cpd_to_vertices_mapper[mesh_vertices[i].cpd_id] = i;
        }
    }
}

void Mesh::create_cpd_out_list(Graph g ,const std::vector<int>& cpd_mapper) {
    std::vector<int> mp = invert_permutation(cpd_mapper);
    g.resort_graph(mp);
    cpd_out_vertices.resize(g.number_of_vertices);
    cpd_distance_cost.resize(g.number_of_vertices);
    for(int i = 0; i < cpd_out_vertices.size();i++){
        for(int j = g.vertices[i]; j < g.vertices[i+1];++j){
            cpd_out_vertices[i].push_back(g.out_vertices[j]);
            cpd_distance_cost[i].push_back(g.distance_cost[j]);
        }
        cpd_out_vertices[i].shrink_to_fit();
        cpd_distance_cost[i].shrink_to_fit();
    }

}


int Mesh::get_map_height(){
    return height;
}

int Mesh::get_map_width(){
    return width;
}

}

