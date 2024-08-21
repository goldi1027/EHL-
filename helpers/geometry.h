#pragma once
#include "point.h"
#include <vector>
#include <algorithm>
/*
 * This file consists of many functions that uses geometry and maths to aid in solving pathfinding problems

 */
namespace polyanya
{

void line_intersect_time(const Point& a, const Point& b,
                         const Point& c, const Point& d,
                         double& ab_num, double& cd_num, double& denom);

// Returns the line intersect between ab and cd as fast as possible.
// Uses a and b for the parameterisation.
// ASSUMES NO COLLINEARITY.
inline Point line_intersect(const Point& a, const Point& b,
                            const Point& c, const Point& d)
{
    const Point ab = b - a;
    return a + ab * (((c - a) * (d - a)) / (ab * (d - c)));
}

inline bool onsegment(const Point& c, const Point& s_1, const Point& s_2){
    double seg_distance = s_1.distance(s_2);
    return s_1.distance(c) < seg_distance && s_2.distance(c) < seg_distance;
}

enum struct ZeroOnePos
{
    LT_ZERO,  // n < 0
    EQ_ZERO,  // n = 0
    IN_RANGE, // 0 < n < 1
    EQ_ONE,   // n = 1
    GT_ONE,   // n > 1
};

enum struct SegIntPos // segment intersect position
{
  DISJOINT,
  INTERSECT,
  OVERLAP
};

// Returns where num / denom is in the range [0, 1].
inline ZeroOnePos line_intersect_bound_check(
    const double num, const double denom
)
{
    // Check num / denom == 0.
    if (std::abs(num) < EPSILON)
    {
        return ZeroOnePos::EQ_ZERO;
    }
    // Check num / denom == 1.
    // Note: Checking whether it is accurately near 1 requires us to check
    // |num - denom| < EPSILON^2 * denom
    // which requires a multiplication. Instead, we can assume that denom
    // is less than 1/EPSILON and check
    // |num - denom| < EPSILON
    // instead. This is less accurate but faster.
    if (std::abs(num - denom) < EPSILON)
    {
        return ZeroOnePos::EQ_ONE;
    }

    // Now we finally check whether it's greater than 1 or less than 0.
    if (denom > 0)
    {
        if (num < 0)
        {
            // strictly less than 0
            return ZeroOnePos::LT_ZERO;
        }
        if (num > denom)
        {
            // strictly greater than 1
            return ZeroOnePos::GT_ONE;
        }
    }
    else
    {
        if (num > 0)
        {
            // strictly less than 0
            return ZeroOnePos::LT_ZERO;
        }
        if (num < denom)
        {
            // strictly greater than 1
            return ZeroOnePos::GT_ONE;
        }
    }
    return ZeroOnePos::IN_RANGE;
}

// Given two points a, b and a number t, compute the point
//  a + (b-a) * t
inline Point get_point_on_line(const Point& a, const Point& b, const double t)
{
    return a + (b - a) * t;
}
Point reflect_point(const Point& p, const Point& l, const Point& r);
inline Point perp_point(const Point& r, const Point& a, const Point& b);

enum struct Orientation
{
	CCW,       // counterclockwise
	COLLINEAR, // collinear
	CW,        // clockwise

};


inline Orientation get_orientation(
    const Point& a, const Point& b, const Point& c
)
{
    const double cross = (b - a) * (c - b);
    if (std::abs(cross) < EPSILON)
    {
        return Orientation::COLLINEAR;
    } else if (cross > 0)
    {
        return Orientation::CCW;
    }
    else
    {
        return Orientation::CW;
    }
}

    enum struct ObstacleLocation
    {
        LEFT,       // counterclockwise
        TOP, // collinear
        RIGHT,        // clockwise

    };


inline ObstacleLocation get_obstacle_side(
        const Point& root, const Point& vertex, const Point& obstacle_vertex1,  const Point& obstacle_vertex2
)
{
    const Orientation&  ov1_ori = get_orientation(root,vertex,obstacle_vertex1);
    const Orientation&  ov2_ori = get_orientation(root,vertex,obstacle_vertex2);
    switch(ov1_ori){
        case Orientation::CCW:
            switch(ov2_ori){
                case Orientation::CCW:
                    return ObstacleLocation::LEFT;
                case Orientation::COLLINEAR:
                    return ObstacleLocation::LEFT;
                case Orientation::CW:
                    return ObstacleLocation::TOP;
            }

        case Orientation::CW:
            switch(ov2_ori){
                case Orientation::CCW:
                    return ObstacleLocation::TOP;
                case Orientation::COLLINEAR:
                    return ObstacleLocation::RIGHT;
                case Orientation::CW:
                    return ObstacleLocation::LEFT;
            }
        case Orientation::COLLINEAR:
            switch(ov2_ori){
                case Orientation::CCW:
                    return ObstacleLocation::LEFT;
                case Orientation::CW:
                    return ObstacleLocation::RIGHT;
                case Orientation::COLLINEAR:
                    assert(false);
                    return ObstacleLocation::RIGHT;
            }
    }
}

    inline void LoadMap(const char *fname, std::vector<std::pair<bool,int>> &map, int &width, int &height)
    {
        FILE *f;
        f = fopen(fname, "r");
        if (f)
        {
            fscanf(f, "type octile\nheight %d\nwidth %d\nmap\n", &height, &width);
            map.resize(height*width);
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    char c;
                    do {
                        fscanf(f, "%c", &c);
                    } while (isspace(c));

                    map[y*width+x].first = (c == '.' || c == 'G' || c == 'S');
                    map[y*width+x].second = y;
                }
            }
            fclose(f);
        }
    }

inline bool is_collinear(const Point& a, const Point& b, const Point& c)
{
    return std::abs((b - a) * (c - b)) < EPSILON;
}

// return [-180, 180] / [0, 360]
inline double get_angle(const Point& p, bool is360=false) {
  double res = std::atan2(p.y, p.x) * 180 / PI;
  if (is360 && res < 0)
    res += 360.0;
  return res;
}

inline double get_relative_angle(const Point& p1, const Point& p2, const Point& p3) {
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

inline double get_point_angle(double p1x, double p1y, double p2x,double p2y)
{
    //Make point1 the origin, make point2 relative to the origin so we do point1 - point1, and point2-point1,
    //since we dont need point1 for the equation to work, the equation works correctly with the origin 0,0.
    double deltaY = p2y - p1y;
    double deltaX = p2x - p1x;

    float angleInDegrees = atan2(deltaY, deltaX) * 180 / 3.141;
    if(angleInDegrees < 0){
    angleInDegrees += 360;
}
return angleInDegrees;
}
inline double get_vector_angle(const Point& p1, const Point& p2) {
    Point p = p2-p1;
    double res = std::atan2(p.y, p.x) * 180 / PI;
    if (res < 0)
            res += 360.0;
    return res;
}

inline bool is_intersect(const Point& p0, const Point& p1, const Point& q0, const Point& q1) {
  // from https://stackoverflow.com/questions/563198/whats-the-most-efficent-way-to-calculate-where-two-line-segments-intersect
    if(q0 == p0 || q0 == p1 || q1==p0 || q1 == p1){
        return true;
    }
  Point r = p1 - p0;
  Point s = q1 - q0;
  if (fabs(r*s) < EPSILON && fabs((q0-p0)*r) > EPSILON) return false;
  if (fabs(r*s) >= EPSILON) {
    double t = (q0-p0)*s / (r*s);
    double u = (p0-q0)*r / (s*r);
    if (0<=t && t<=1 && 0<=u && u<=1) return true;
    return false;
  }
  return false;
}








inline bool same_side(const Point& v1, const Point& v2, const Point& l, const Point& r) {
  const Point lr = r - l;
  const Point lv1= v1 - l;
  const Point lv2 = v2 - l;
  if ((lv1 * lr > 0) == (lv2 * lr > 0)) 
    return true;
  return false;
}

inline bool onSegment(const Point& p, const Point& q, const Point& r)
    {
        if (q.x <= fmax(p.x, r.x) && q.x >= fmin(p.x, r.x) &&
            q.y <= fmax(p.y, r.y) && q.y >= fmin(p.y, r.y))
            return true;

        return false;
    }

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
inline int orientation(const Point& p, const Point& q, const Point& r)
    {
        // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
        // for details of below formula.
        int val = (q.y - p.y) * (r.x - q.x) -
                  (q.x - p.x) * (r.y - q.y);

        if (val == 0) return 0;  // colinear

        return (val > 0)? 1: 2; // clock or counterclock wise
    }

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
 inline bool doIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2)
    {
        // Find the four orientations needed for general and
        // special cases
        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);

        // General case
        if (o1 != o2 && o3 != o4)
            return true;

        // Special Cases
        // p1, q1 and p2 are colinear and p2 lies on segment p1q1
        if (o1 == 0 && onSegment(p1, p2, q1)) return true;

        // p1, q1 and q2 are colinear and q2 lies on segment p1q1
        if (o2 == 0 && onSegment(p1, q2, q1)) return true;

        // p2, q2 and p1 are colinear and p1 lies on segment p2q2
        if (o3 == 0 && onSegment(p2, p1, q2)) return true;

        // p2, q2 and q1 are colinear and q1 lies on segment p2q2
        if (o4 == 0 && onSegment(p2, q1, q2)) return true;

        return false; // Doesn't fall in any of the above cases
    }


    inline Point lineLineIntersection(Point A, Point B, Point C, Point D)
    {
        // Line AB represented as a1x + b1y = c1
        double a1 = B.y - A.y;
        double b1 = A.x - B.x;
        double c1 = a1*(A.x) + b1*(A.y);

        // Line CD represented as a2x + b2y = c2
        double a2 = D.y - C.y;
        double b2 = C.x - D.x;
        double c2 = a2*(C.x)+ b2*(C.y);

        double determinant = a1*b2 - a2*b1;

        if (determinant == 0)
        {
            // The lines are parallel. This is simplified
            // by returning a pair of FLT_MAX
            return Point{-1,-1};
        }
        else
        {
            double x = (b2*c1 - b1*c2)/determinant;
            double y = (a1*c2 - a2*c1)/determinant;
            return Point{x, y};
        }
    }

    inline bool raytTolineIntersection(Point p1, Point p2, Point q1, Point q2){
        Point intersection = lineLineIntersection(p1, p2, q1, q2);
        if(intersection == Point {-1,-1}){
            //not intersect
            if(get_orientation(q1,q2,p1) == Orientation::COLLINEAR){
                return onSegment(q1,p1,q2);
            }
            if(get_orientation(q1,q2,p2) == Orientation::COLLINEAR){
                return onSegment(q1,p2,q2);
            }
            return false;
        }else{
            if(onSegment(q1,intersection,q2)){
                if(onSegment(p1,intersection,p2)){
                    return true;
                }else if (intersection.distance(p2) <= intersection.distance(p1) + EPSILON){
                    return true;
                }
            }

        }
        return false;


    }
    int inSegment(const Point& p, const Point& s1, const Point& s2);
SegIntPos intersect2D_2Segments(const Point& p0, const Point& p1, const Point& q0, const Point& q1, Point& I0, Point& I1);
SegIntPos intersect2D_2Segments(const Point& p0, const Point& p1, const Point& q0, const Point& q1, Point& I0, Point& I1);
bool is_permutation(const std::vector<int>&p);
std::vector<int> invert_permutation(const std::vector<int>&p);
}
