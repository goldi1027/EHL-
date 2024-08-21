/*
 * This file consists of many functions that uses geometry and maths to aid in solving pathfinding problems

 */

#include "geometry.h"
#include "point.h"
#include "consts.h"
#include <cmath>
#include <cassert>

namespace polyanya
{

// Given two line segments ab and cd defined as:
//   ab(t) = a + (b-a) * t
//   cd(t) = c + (d-c) * t
// find ab_num, cd_num, denom such that
//   ab(ab_num / denom) = cd(cd_num / denom).
// If denom is close to 0, it will be set to 0.
void line_intersect_time(const Point& a, const Point& b,
                         const Point& c, const Point& d,
                         double& ab_num, double& cd_num, double& denom)
{
    denom = (b - a) * (d - c);
    if (std::abs(denom) < EPSILON)
    {
        denom = 0.0; // to make comparison easy
        ab_num = 1;
        cd_num = 1;
    }
    else
    {
        const Point ac = c - a;
        ab_num = ac * (d - a);
        cd_num = ac * (b - c);

        #ifndef NDEBUG
        // If we're debugging, double check that our results are right.
        const Point ab_point = get_point_on_line(a, b, ab_num / denom);
        const Point cd_point = get_point_on_line(c, d, cd_num / denom);
        assert(ab_point == cd_point);
        #endif
    }
}

// Reflects the point across the line lr.
Point reflect_point(const Point& p, const Point& l, const Point& r)
{
    // I have no idea how the below works.
    const double denom = r.distance_sq(l);
    if (std::abs(denom) < EPSILON)
    {
        // A trivial reflection.
        // Should be p + 2 * (l - p) = 2*l - p.
        return 2 * l - p;
    }
    const double numer = (r - p) * (l - p);

    // The vector r - l rotated 90 degrees counterclockwise.
    // Can imagine "multiplying" the vector by the imaginary constant.
    const Point delta_rotated = {l.y - r.y, r.x - l.x};
    return p + (2.0 * numer / denom) * delta_rotated;
}

inline Point perp_point(const Point& r, const Point& a, const Point& b) {
  if (std::abs(a.x - b.x) <= EPSILON && std::abs(a.y - b.y) <= EPSILON) {
    return Point{a.x, a.y};
  } else {
    Point p = b - a;
    double u = p.dot(r - a) / p.normal2();
    Point res = r + u * p;
    return res;
  }
}

// intersect2D_2Segments(): find the 2D intersection of 2 finite segments
//    Input:  two finite segments S1(p0, p1) and S2(q0, q1)
//    Output: I0 = intersect point (when it exists)
//            I1 =  endpoint of intersect segment [I0,I1] (when it exists)
//    Return: disjoint (no intersect)
//            intersect  in unique point I0
//            overlap  in segment from I0 to I1
// from http://geomalgorithms.com/a05-_intersect-1.html#Line-Intersections
SegIntPos intersect2D_2Segments(const Point& p0, const Point& p1, const Point& q0, const Point& q1,
                          Point& I0, Point& I1) {
  Point u = p1 - p0;
  Point v = q1 - q0;
  Point w = p0 - q0;
  double D = u * v;
  if (fabs(D) < EPSILON) { // S1 and S2 are parallel
    if (fabs(u * w) > EPSILON || fabs(v * w) > EPSILON) { // they are NOT collinear
      return SegIntPos::DISJOINT;
    }
    double du = u.dot(u);
    double dv = v.dot(v);
    if (fabs(du) < EPSILON && fabs(dv) < EPSILON) { // both are points
      if (p0.distance(q0) > EPSILON) // distinct points
        return SegIntPos::DISJOINT;
      I0 = Point{p0.x, p0.y};
      return SegIntPos::INTERSECT;
    }
    if (fabs(du) < EPSILON) { // S1 is a single point
      if (inSegment(p0, q0, q1) == 0) // not in S2
        return SegIntPos::DISJOINT;
      I0 = Point{p0.x, p0.y};
      return SegIntPos::INTERSECT;
    }
    if (fabs(dv) < EPSILON) { // S2 is a single point
      if (inSegment(q0, p0, p1) == 0) // not in S1
        return SegIntPos::DISJOINT;
      I0 = Point{q0.x, q0.y};
      return SegIntPos::INTERSECT;
    }
    // they are collinear segments - get overlap or not
    double t0, t1;
    Point w2 = p1 - q0;
    if (fabs(v.x) > EPSILON) { // end of s1 in eqn for s2
      t0 = w.x / v.x;
      t1 = w2.x / v.x;
    } else {
      t0 = w.y / v.y;
      t1 = w2.y / v.y;
    }
    if (t0 > t1) {
      std::swap(t0, t1);
    }
    if (t0 > 1.0 + EPSILON || t1 < -EPSILON) {
      return SegIntPos::DISJOINT; // NOT overlap
    }
    t0 = t0<0?0: t0;
    t1 = t1>1? 1: t1;
    if (fabs(t0 - t1) < EPSILON) { // the intersect is a point
      I0 = q0 + t0 * v;
      return SegIntPos::INTERSECT;
    }

    // they overlap in a valid subsegment
    I0 = q0 + t0 * v;
    I1 = q0 + t1 * v;
    return SegIntPos::OVERLAP;
  }

  // segments are skew and may intersect in a point
  // get the intersect parameter for s1
  double sI = (v * w) / D;
  if (sI <= -EPSILON || sI >= 1 + EPSILON) // not intersect with s1
    return SegIntPos::DISJOINT;

  // get the intersect parameter for s2
  double tI = u * w / D;
  if (tI <= -EPSILON || tI >= 1 + EPSILON) // not intersect with s2
    return SegIntPos::DISJOINT;

  I0 = p0 + sI * u;
  return SegIntPos::INTERSECT;
}

// inSegment(): determine if a point is inside a segment
//    Input:  a point P, and a collinear segment S
//    Return: 1 = P is inside S
//            0 = P is  not inside S
int inSegment( const Point& P, const Point& s0, const Point& s1) {
    if (s0.x != s1.x) {    // S is not  vertical
        if (s0.x <= P.x && P.x <= s1.x)
            return 1;
        if (s0.x >= P.x && P.x >= s1.x)
            return 1;
    }
    else {    // S is vertical, so test y  coordinate
        if (s0.y <= P.y && P.y <= s1.y)
            return 1;
        if (s0.y >= P.y && P.y >= s1.y)
            return 1;
    }
    return 0;
}
    bool is_permutation(const std::vector<int>&p){
        std::vector<bool>found(p.size(), false);
        for(unsigned x:p){
            if(x >= p.size())
                return false;
            if(found[x])
                return false;
            found[x] = true;
        }
        return true;
    }
    std::vector<int> invert_permutation(const std::vector<int>&p){
        assert(is_permutation(p) && "p must be a permutation");

        std::vector<int> inv_p(p.size());
        for(int i=0; i<p.size(); ++i)
            inv_p[p[i]] = i;

        return inv_p; // NVRO
    }

}
