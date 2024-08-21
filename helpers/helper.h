#ifndef LEARNED_HEURISTIC_HELPER_H
#define LEARNED_HEURISTIC_HELPER_H
#include <cmath>

bool isInside(int x1, int y1, int x2, int y2, int x3, int y3, int x, int y)
{
    /* Calculate area of triangle ABC */
    float A = abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2.0);

    /* Calculate area of triangle PBC */
    float A1 = abs((x*(y2-y3) + x2*(y3-y)+ x3*(y-y2))/2.0);
    //area (x, y, x2, y2, x3, y3);

    /* Calculate area of triangle PAC */
    float A2 = abs((x1*(y-y3) + x*(y3-y1)+ x3*(y1-y))/2.0);
    //area (x1, y1, x, y, x3, y3);

    /* Calculate area of triangle PAB */
    float A3 = abs((x1*(y2-y) + x2*(y-y1)+ x*(y1-y2))/2.0);
    //area (x1, y1, x2, y2, x, y);

    /* Check if sum of A1, A2 and A3 is same as A */
    return (A == A1 + A2 + A3);
}

float get_angles(double p1x, double p1y, double p2x,double p2y)
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


double get_angles_between_points(double p1x, double p1y, double p2x,double p2y, double p3x,double p3y){
    return atan2(p3y- p1y, p3x - p1x) - atan2(p2y- p1y, p2x- p1x);
}

#endif //LEARNED_HEURISTIC_HELPER_H
