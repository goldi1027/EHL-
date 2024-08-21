/*
An Experimental Study on Hub Labeling based Shortest Path Algorithms [Experiments and Analyses]

Authors: Ye Li, Leong Hou U, Man Lung Yiu, Ngai Meng Kou
Contact: yb47438@umac.mo
Affiliation: University of Macau

The MIT License (MIT)

Copyright (c) 2016 University of Macau

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
#pragma once

#ifndef TIME_UTIL_H
#define TIME_UTIL_H

#include <sys/time.h>
#include<map>
#include <cmath>

namespace time_util {

	//double GetCurrentTimeSec() { return 0; }

	
	double GetCurrentTimeSec() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }

	void save_dis_stat(string dis_stat_file, vector<pair<int, int> > dis_meet, vector<pair<int, int> > dis_meet_dis1, vector<pair<int, int> > dis_meet_dis2){
	
		ofstream ofs(dis_stat_file);
		
		for(int i = 0; i < dis_meet.size(); ++i)
			ofs << dis_meet[i].first << "\t" << dis_meet[i].second << "\t" << dis_meet_dis1[i].second << "\t" << dis_meet_dis2[i].second << endl;		
				
		ofs.close();	
	}

	void save_stat(string stat_file, double load_time, double ordering_time, double labeling_time, double avg_size, vector<double>& each_iteration, vector<double>& pruning_power) {
		ofstream ofs(stat_file);

		ofs << numOfVertices << "\t" << numOfEdges		<< endl;
		ofs << load_time << endl;
		ofs << ordering_time << endl;
		ofs << labeling_time << endl;
		ofs << avg_size << endl;

		for (NodeID v = 0; v < each_iteration.size(); ++v) {
		ofs << each_iteration[v] << endl;
		}

		ofs << "pruning_power" << endl;
		for (NodeID v = 0; v < pruning_power.size(); ++v) {
		ofs << pruning_power[v] << endl;
		}
		

		ofs.close();
	}
	
	void save_stat_betweenness(string stat_file, double init_time, double selecting_time, double labeling_time, double updating_time, double adding_time, double avg_size, vector<double>& each_iteration, vector<NodeID>& pruning_power) {
		ofstream ofs(stat_file);
		double total = init_time + selecting_time + labeling_time + updating_time + adding_time;
		ofs << numOfVertices << "\t" << numOfEdges << endl;
		ofs << init_time/total << endl;
		ofs << selecting_time/total << endl;
		ofs << labeling_time/total << endl;
		ofs << updating_time/total << endl;
		ofs << adding_time/total << endl;
		ofs << total << endl;
		ofs << avg_size << endl;

		for (NodeID v = 0; v < each_iteration.size(); ++v) {
		ofs << each_iteration[v] << endl;
		}

		ofs << "pruning_power" << endl;
		for (NodeID v = 0; v < pruning_power.size(); ++v) {
		ofs << pruning_power[v] << endl;
		}
		

		ofs.close();
	}

	
}
#endif