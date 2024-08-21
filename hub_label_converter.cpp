/*
 * Converts original hub label file format into usable format for computing EHL
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iomanip>


using namespace std;
void load_data(string dimacs_file,vector<unsigned>& head,vector<unsigned>& tail,vector<double>& weight){

    cout << "Loading data ... " << flush;

    ifstream in(dimacs_file);
    if(!in)
        throw runtime_error("Can not open \""+dimacs_file+"\"");

    string line;
    unsigned line_num = 0;
    unsigned next_arc = 0;

    unsigned node_count, arc_count;

    bool was_node_read = false;
    bool was_arc_read = false;
    int i = 0;
    while(std::getline(in, line)){
        ++line_num;

        std::istringstream lin(line);
        if(i == 0){
            //std::string p, sp;
            if(!(lin >> node_count))
                throw std::runtime_error("Can not parse node count in file.");}

        if(i == 1){
            if (!(lin >> arc_count))
                throw std::runtime_error("Can not parse arc count in file.");

            tail.resize(arc_count);
            head.resize(arc_count);
            weight.resize(arc_count);
        }
        if(i>1){
            int h, t;
            double w;
            if(!(lin >> t >> h >> w))
                throw std::runtime_error("Can not parse line num "+std::to_string(line_num)+" \""+line+"\" in dimacs file.");





            if(next_arc < arc_count){
                head[next_arc] = h;
                tail[next_arc] = t;
                weight[next_arc] = w;
            }
            ++next_arc;
        }
        ++i;
    }

    if(next_arc != arc_count)
        throw std::runtime_error("The arc count in the header ("+to_string(arc_count)+") does not correspond with the actual number of arcs ("+to_string(next_arc)+").");

    cout << "done" << endl;


}



int main(int argc, char*argv[]){

    try{
        string dimacs_distance_file;
        string txt_file;
//        string weight_file;



        dimacs_distance_file = argv[1];
        txt_file = argv[2];


        vector<unsigned> distance_head,distance_tail;
        vector<double>distance;

        load_data(dimacs_distance_file,distance_head,distance_tail,distance);

        std::ofstream myFile(txt_file);
        bool should_remove = true;
        for(int i = 0 ; i < distance_head.size(); i ++){
            myFile<<std::fixed<<setprecision(8)<<distance_head[i]<<" "<<distance_tail[i]<<" "<<distance[i]<<"\n";
        }
        myFile.close();


        cout << "done" << endl;


    }catch(exception&err){
        cerr << "Stopped on exception : " << err.what() << endl;
    }
}
