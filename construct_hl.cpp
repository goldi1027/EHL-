/*
 * Computes hub label based on a given visibility graph
 */


#include "graph.h"
#include "coverage_ordering_path.h"
#include <stdexcept>
#include <iomanip>
#include "dijkstra.h"
#include <iostream>
#include <searchinstance.h>

using namespace std;
namespace pl = polyanya;
pl::MeshPtr mp;
pl::SearchInstance* si;

void build_shp_undirect_graph(const char *graphFileName,const char * labelFileName,const char * output, const char *file_name, string directory_name){
    WGraph wgraph;
    wgraph.load_graph(graphFileName);
    double avg_degree;
    double max_degree = 0;
    Coverage_Ordering_Path coverage_ordering(wgraph, avg_degree, max_degree);

    string orderFileName(labelFileName);
    orderFileName.append(".order");
    coverage_ordering.save_rank(orderFileName.c_str());

    string labelFile(labelFileName);
    labelFile.append(".label");
    //save hub label file
    coverage_ordering.plabels.save_labels(labelFile.c_str());

}


int main(int argc, char *argv[]) {

    build_shp_undirect_graph(argv[1],argv[2],argv[3], argv[4], argv[5]);

    return 0;
}