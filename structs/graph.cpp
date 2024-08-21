#include "graph.h"
#include <stack>
#include <iostream>
namespace polyanya {

    void Graph::load_graph(const string &file_name) {
        ifstream in(file_name);
        vector<vector<double>> lines;
        number_of_vertices = 0 ;
        if (in) {
            string line;
            int i = 0;
            while (getline(in, line)) {
                if(i == 0){
                    number_of_vertices = stoi(line);
                }
                if(i >1) {
                    stringstream sep(line);
                    string field;
                    lines.push_back(vector<double>());
                    while (getline(sep, field, ' ')) {
                        lines.back().push_back(stod(field));
                    }
                }
                i++;
            }
        }

        vector<vector<int>> visibility_graph ;
        vector<vector<double>> graph_weight;

        visibility_graph.resize(number_of_vertices);
        graph_weight.resize(number_of_vertices);
        number_of_edges = 0;
        for( auto line : lines){
            unsigned head = line[0];
            unsigned tail = line[1];
            double weight = line[2];
            visibility_graph[head].push_back(tail);
            graph_weight[head].push_back(weight);
            number_of_edges++;
        }
        vertices = vector<int>(number_of_vertices+1);
        out_vertices = vector<int>(number_of_edges);
        distance_cost = vector<double>(number_of_edges);
        int start_index = 0;
        int max_degree =0;
        for(int i = 0; i <visibility_graph.size(); i++){
            vertices[i] = start_index;
            if(max_degree < visibility_graph[i].size()){
                max_degree = visibility_graph[i].size();
            }
            for(int j = 0; j < visibility_graph[i].size();j++){
                out_vertices[start_index + j] = visibility_graph[i][j];
                distance_cost[start_index + j] = graph_weight[i][j];

            }
            start_index += visibility_graph[i].size();
        }
        vertices[number_of_vertices] = start_index;
        if(max_degree > 4096){
            cout<<"Degree is larger than CPD setting! Probably you can resolve this by cutting necessary first move during the preprocessing "<<endl;
        }


    }


// may use this to generate query:
    bool Graph::is_connected(int s, int t) {
        vector<bool> visited(number_of_vertices, false);

        // Create a stack for DFS
        stack<int> stack;

        // Push the current source node.
        stack.push(s);

        while (!stack.empty()) {
            // Pop a vertex from stack and print it
            s = stack.top();
            stack.pop();

            // Stack may contain same vertex twice. So
            // we need to print the popped item only
            // if it is not visited.
            if (!visited[s]) {
//                cout << s << " ";
                visited[s] = true;
                if (s == t) {
                    return true;
                }
            }

            // Get all adjacent vertices of the popped vertex s
            // If a adjacent has not been visited, then push it
            // to the stack.
            for (int a = vertices[s]; a < vertices[s + 1]; ++a) {
                if (!visited[out_vertices[a]]) {
                    stack.push(out_vertices[a]);
                }

            }
        }
        return false;
    }


    vector<int> Graph::generate_DFS_ordering() {

        vector<int> DFS_ordering(number_of_vertices,0);
        // Mark all the vertices as not visited

        vector<bool> was_pushed(number_of_vertices, false);
        vector<bool> was_popped(number_of_vertices, false);
        std::vector<int>next_out = vertices;
        stack<int>to_visit;
        int next_id = 0;
        for(int source_node=0; source_node<number_of_vertices; ++source_node){
            if(!was_pushed[source_node]){

                to_visit.push(source_node);
                was_pushed[source_node] = true;

                while(!to_visit.empty()){
                    int x = to_visit.top();
                    to_visit.pop();
                    if(!was_popped[x]){
                        DFS_ordering[x]= next_id++;
                        was_popped[x] = true;
                    }

                    while(next_out[x] != vertices[x+1]){
                        int y = out_vertices[next_out[x]];
                        if(was_pushed[y])
                            ++next_out[x];
                        else{
                            was_pushed[y] = true;
                            to_visit.push(x);
                            to_visit.push(y);
                            ++next_out[x];
                            break;
                        }
                    }
                }
            }
        }
        DFS_ordering = invert_permutation( DFS_ordering);
        assert(DFS_ordering.size() == number_of_vertices);
        return DFS_ordering;
    }


    void Graph::resort_graph(const vector<int> &ordering) {
        //resort the graph based on certain ordering
        //This function should be used for building CPD only.
        //only resorted vertices, out_vertices, distance_cost.

        //inverted mapper, map current v_id to ordering.
        vector<int> inv_ordering = invert_permutation(ordering);

        vector<int> tmp_vertices;
        vector<int> tmp_out_vertices;
        vector<double> tmp_distance_cost;
        assert(ordering.size() == number_of_vertices);
        int start_index = 0;
        for (int i = 0; i < number_of_vertices; i++) {
            tmp_vertices.push_back(start_index);
            int v_id = ordering[i];
            for (int arc = vertices[v_id]; arc < vertices[v_id + 1]; ++arc) {
                tmp_out_vertices.push_back(inv_ordering[out_vertices[arc]]);
                tmp_distance_cost.push_back( distance_cost[arc]);
                start_index++;
            }
        }
        tmp_vertices.push_back(start_index);

        assert(tmp_vertices.size() == vertices.size());
        assert(tmp_out_vertices.size() == out_vertices.size());
        assert(tmp_distance_cost.size() == distance_cost.size());

        vertices = tmp_vertices;
        out_vertices = tmp_out_vertices;
        distance_cost = tmp_distance_cost;
    }

}