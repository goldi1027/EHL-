/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <id_queue.h>
#include <constants.h>
#include <timestamp_flag.h>
#include <vector>
#include <iostream>
#include <set>
using namespace std;
namespace polyanya{
class Dijkstra{
public:

	Dijkstra(const vector<vector<unsigned>>&visibility_graph, const vector<vector<double>>&graph_weight):
		queue(visibility_graph.size()),
        was_distance_set(visibility_graph.size()),
        visibility_graph(&visibility_graph),
        graph_weight(&graph_weight){
		assert(!visibility_graph.empty());
        dist.resize(visibility_graph.size());
	}

	Dijkstra&reset(){
		queue.clear();
		return *this;
	}

    void reach(const unsigned & next_vertex, double d){
        if(was_distance_set.is_set(next_vertex)){
            if(d <  dist[next_vertex]){
                queue.decrease_key({next_vertex, d});
                dist[next_vertex] = d;
            }
        }else{
            queue.push({next_vertex,d});
            dist[next_vertex] = d;
            was_distance_set.set(next_vertex);
        }

    }
    const std::vector<double>& get_distance_to_all(int source_node){
	    std::fill(dist.begin(),dist.end(),-1);
        queue.clear();
        was_distance_set.reset_all();
        dist[source_node] = 0;
        was_distance_set.set(source_node);

        for(int i = 0; i<(*visibility_graph)[source_node].size();i++){
            reach((*visibility_graph)[source_node][i], (*graph_weight)[source_node][i]);
        }

        while(!queue.empty()){
            auto x = queue.pop();
            for(int i = 0; i<(*visibility_graph)[x.id].size();i++ ){
                double weight = (*graph_weight)[x.id][i];
                reach((*visibility_graph)[x.id][i], x.key + weight);
            }
        }

        return dist;
    }
private:
    std::vector<double>dist;
	MinIDQueue queue;
    TimestampFlags was_distance_set;
    const vector<vector<unsigned>>*visibility_graph;
    const vector<vector<double>>*graph_weight;


};


}

#endif
