#pragma once
#include <vector>

namespace polyanya
{

    struct Edge
    {
//        int id = 0 ;
        // "int" here means a pair index.
        std::pair<int,int> vertices;
        std::pair<int,int> polygons;
    };
}