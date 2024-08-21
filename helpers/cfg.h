/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 */
#pragma once

#include "getopt.h"

#include <iostream>
#include <map>

namespace warthog
{

namespace util
{

typedef struct option param;

class cfg
{
    public:
        cfg();
        ~cfg();

        void 
        parse_args(int argc, char** argv, warthog::util::param params[]);

        std::string
        get_param_value(std::string);


        void
        print(std::ostream& out);

    private:
        // no copy
        cfg(const cfg& other) {}
        cfg& operator=(const cfg& other) { return *this; }

        std::map<std::string, std::string> params_;
};

}
}
