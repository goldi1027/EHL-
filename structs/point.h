/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#pragma once
#include "consts.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "polygon.h"
#include <unordered_set>

namespace polyanya
{

// An (x, y) pair.
    struct Point
    {
        double x, y;
        int vertexId = -1;

        bool operator<(const Point& other) const {
          return x <= other.x && y <= other.y;
        }

        bool operator==(const Point& other) const
        {
            return (std::abs(x - other.x) < EPSILON) &&
                   (std::abs(y - other.y) < EPSILON);
        }

        bool operator!=(const Point& other) const
        {
            return !((*this) == other);
        }

        Point operator+(const Point& other) const
        {
            return {x + other.x, y + other.y};
        }

        Point operator-() const
        {
            return {-x, -y};
        }

        Point operator-(const Point& other) const
        {
            return {x - other.x, y - other.y};
        }

        // Cross product.
        // Returns the z component (as we are working in 2D).
        double operator*(const Point& other) const
        {
            return x * other.y - y * other.x;
        }

        Point operator*(const double& mult) const
        {
            return {mult * x, mult * y};
        }

        friend Point operator*(const double& mult, const Point& p)
        {
            return p * mult;
        }

        friend std::ostream& operator<<(std::ostream& stream, const Point& p)
        {
            return stream << "(" << p.x << ", " << p.y << ")";
        }

        double distance_sq(const Point& other) const
        {
            #define square(x) (x)*(x)
            return square(x-other.x) + square(y-other.y);
            #undef square
        }

        //euclidean distance?
        double distance(const Point& other) const
        {
            return std::sqrt(this->distance_sq(other));
        }

        double dot(const Point& other) const {
          return this->x * other.x + this->y * other.y;
        }

        double normal() {
          return std::sqrt(this->x * this->x + this->y * this->y);
        }

        double normal2() {
          return this->x * this->x + this->y * this->y;
        }

        double distance_to_seg(const Point& l, const Point& r) const {
          Point b = r - l, p = *this - l;
          if (std::abs(b.x) < EPSILON && std::abs(b.y) < EPSILON)
            return this->distance(l);
          double t = b.dot(p) / b.normal2();
          Point d;
          if (t < 0) {
            d = (*this) - l;
          }
          else if (t > 1) {
            d = (*this) - r;
          }
          else {
            d = (*this) - (l + t * b);
          }
          return d.normal();
        }

        double min_distance_to_grid(Point min, Point max){
            double rx = (max.x + min.x)/2;
            double ry = (max.y + min.y)/2;
            double width = max.x-min.x;
            double height = max.y-min.y;

            double v1 = fabs(this->x  - rx) - width / 2;
            double v2 = fabs(this->y  - ry) - height / 2;

            double dx = (v1>0 ?  v1 : 0) ;
            double dy = (v2>0 ?  v2 : 0) ;
            return std::sqrt(dx * dx + dy * dy);
        }

        double max_distance_to_grid(Point min, Point max){
            double far_x =  fabs(this->x  - min.x) >=fabs(this->x  - max.x) ? min.x :max.x;
            double far_y =  fabs(this->y  - min.y) >=fabs(this->y  - max.y) ? min.y :max.y;
            return std::sqrt((far_x-this->x) *(far_x-this->x) + (far_y-this->y)*(far_y-this->y));
        }

        double distance_to_grid(int x, int y) const{
            double min_x = x - 0.5;
            double max_x = x + 0.5;
            double min_y = y - 0.5;
            double max_y = y + 0.5;

            double rx = (max_x + min_x)/2;
            double ry = (max_y + min_y)/2;
            double width = max_x-min_x;
            double height = max_y-min_y;

            double v1 = fabs(this->x  - rx) - width / 2;
            double v2 = fabs(this->y  - ry) - height / 2;

            double dx = (v1>0 ?  v1 : 0) ;
            double dy = (v2>0 ?  v2 : 0) ;
            return std::sqrt(dx * dx + dy * dy);
        }

        double max_distance_to_grid(int x, int y) const{
            double min_x = x - 0.5;
            double max_x = x + 0.5;
            double min_y = y - 0.5;
            double max_y = y + 0.5;

            double far_x =  fabs(this->x  - min_x) >=fabs(this->x  - max_x) ? min_x :max_x;
            double far_y =  fabs(this->y  - min_y) >=fabs(this->y - max_y) ? min_y :max_y;

            return std::sqrt((far_x-this->x) *(far_x-this->x) + (far_y-this->y)*(far_y-this->y));

            double rx = (max_x + min_x)/2;
            double ry = (max_y + min_y)/2;
            double width = max_x-min_x;
            double height = max_y-min_y;

            double v1 = fabs(this->x  - rx) - width / 2;
            double v2 = fabs(this->y  - ry) - height / 2;

            double dx = (v1>0 ?  v1 : 0) ;
            double dy = (v2>0 ?  v2 : 0) ;
            return std::sqrt(dx * dx + dy * dy);
        }

        double distance_to_rectangle(Polygon p){
            double min_x = p.min_x;
            double max_x = p.max_x;
            double min_y = p.min_y;
            double max_y = p.max_y;

            double rx = (max_x + min_x)/2;
            double ry = (max_y + min_y)/2;
            double width = max_x-min_x;
            double height = max_y-min_y;

            double v1 = fabs(this->x  - rx) - width / 2;
            double v2 = fabs(this->y  - ry) - height / 2;

            double dx = (v1>0 ?  v1 : 0) ;
            double dy = (v2>0 ?  v2 : 0) ;
            return std::sqrt(dx * dx + dy * dy);
        }
    };

}
