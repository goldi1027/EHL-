/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#include "scenario.h"
#include "point.h"
#include <iostream>
#include <vector>
#include <string>

namespace polyanya
{

bool load_scenarios(std::istream& infile, std::vector<Scenario>& out)
{
    std::string version_str;
    double version;
    if (!(infile >> version_str >> version))
    {
        std::cerr << "Couldn't find scenario version." << std::endl;
        return false;
    }
    if (version_str != "version" || version != 1.0)
    {
        std::cerr << "Invalid scenario version." << std::endl;
        return false;
    }
    Scenario s;
    std::string map;
    while (infile >> s.bucket >> map >> s.xsize >> s.ysize >> s.start.x >>
           s.start.y >> s.goal.x >> s.goal.y >> s.gridcost)
    {
        out.push_back(s);
    }
    return true;
}

    bool load_workload_scenarios(std::istream& infile, std::vector<Scenario>& out)
    {
        std::cout << "Loading workload scenario" << std::endl;
        std::string version_str;
        double version;
        if(!infile.good()){
            std::cout << "cannot read file" << std::endl;
        }
        if (!(infile >> version_str >> version))
        {
            std::cerr << "Couldn't find scenario version." << std::endl;
            return false;
        }
        if (version_str != "version" || version != 1.0)
        {
            std::cerr << "Invalid scenario version." << std::endl;
            return false;
        }
        // Skip the first 3000 lines after the header
        int skipLines = 500001;
        std::string tempLine;
        while (skipLines-- > 0 && std::getline(infile, tempLine)) {
            // Just iterating through the lines without processing
        }

        // Read remaining lines into the scenarios
        Scenario s;
        std::string map;
        while (infile >> s.bucket >> map >> s.xsize >> s.ysize >> s.start.x >>
                      s.start.y >> s.goal.x >> s.goal.y >> s.gridcost) {
            out.push_back(s);
        }

        return true;
    }

}
