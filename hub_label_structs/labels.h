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

#ifndef LABELS_H
#define LABELS_H

#include <filesystem>
#include <limits>
#include <climits>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include "graph.h"
#include "paras.h"
#include <stdlib.h>
#include <xmmintrin.h>
//typedef unsigned __int64 BPSeed;
#include <omp.h>
#include<bitset>
#include <cmath>
#include "point.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define INF_WEIGHT SP_Constants::INF_WEIGHT

#ifdef __APPLE__
#include <stdlib.h>
    // wrapper between posix_memalgin and memalign
    void *memalign(size_t blocksize, size_t bytes)
    {
        void *m;
        errno = posix_memalign(&m, blocksize, bytes);
        void* returns = errno ? NULL : m;
        if(returns == NULL){
            std::cout<<"can no allocate memory;"<<std::endl;
        }

        return errno ? NULL : m;
    }
#else
#include <malloc.h>
#endif

struct index_t {
    vector<NodeID> spt_v;
    vector<EdgeWeight> spt_d;

    NodeID size() {
        return spt_v.size();
    }

};

struct index_t_p {
    NodeID* spt_v;
    EdgeWeight* spt_d;
}__attribute__((aligned(64)));  // Aligned for cache lines;


struct two_index_t_p {
    NodeID* spt_v;
    EdgeWeight* spt_d;
    uint8_t* spt_lv;
    EdgeWeight* spt_ld;
}__attribute__((aligned(64)));  // Aligned for cache lines;

struct index_t_path {
    vector<NodeID> spt_v;
    vector<NodeID> spt_p;//parent nodes
    vector<EdgeWeight> spt_d;

    NodeID size() {
        return spt_v.size();
    }

};

struct index_t_path_p {
    NodeID* spt_v;
    NodeID* spt_p;
    EdgeWeight* spt_d;

};

struct query_info {
    NodeID meet_node;
    NodeID search_len;
    double time_cost;
    EdgeWeight distance;
};

template<int kNumBitParallelRoots = 50>
struct index_t_bp {
    NodeID* spt_v;
    EdgeWeight* spt_d;
    EdgeWeight bpspt_d[kNumBitParallelRoots];
    uint64_t bpspt_s[kNumBitParallelRoots][2];
}__attribute__((aligned(64)));  // Aligned for cache lines;


struct token_t {
    NodeID* sptc_v; // sptc_v[0] is the root
    EdgeWeight* sptc_d;	 // |*| = k + 1, sptc_d[0] is the number of children - k
    unsigned char* sptc_fbv; // first-level bit vector
    unsigned char* sptc_sbv; // second-level bit vector
    NodeID* sptc_pathv; // intermediate point for a path
}__attribute__((aligned(64)));

struct raw_label{
    int hub_id;
    double label_distance;
    int label_predecessor;
    int label_first_node;
};


class PLabel {

public:
    vector<index_t_path> index_;
    index_t_path_p* index_p;
    vector<vector<raw_label>> raw_label_list;

    double GetCurrentTimeSec() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }


    PLabel() {
        index_.resize(numOfVertices);
    }

    ~PLabel() {
        Free();
    }

    EdgeWeight query_p(NodeID s, NodeID t) {

        EdgeWeight distance = INF_WEIGHT;
        NodeID meet;

        const index_t_path_p &idx_s = index_p[s];
        const index_t_path_p &idx_t = index_p[t];

        _mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

        for (int i = 0, j = 0; ; ) {
            NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

            if (v1 == numOfVertices) break;  // Sentinel

            if (v1 == v2) {
                EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
                if (td < distance) {
                    distance = td;
                }
                ++i;
                ++j;
            }
            else {
                i += v1 < v2 ? 1 : 0;
                j += v1 > v2 ? 1 : 0;
            }
        }

        return distance;
    }

    void get_label(NodeID s, vector<int> &label_list, vector<double> &label_distance){
        const index_t_path_p &idx = index_p[s];
        _mm_prefetch(&idx.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx.spt_d[0], _MM_HINT_T0);
        for(int i = 0; i < numOfVertices; ++i){
            //if(i>0 and idx.spt_v[i]==NULL){
            //boundary case to stop
            if(idx.spt_v[i] == numOfVertices){
//                label_list.push_back(idx.spt_v[i]);
//                label_distance.push_back(idx.spt_d[i]);
                break;
            }
            //if(idx.spt_d[i] != INF_WEIGHT){
            else{
                label_list.push_back(idx.spt_v[i]);
                label_distance.push_back(idx.spt_d[i]);
                //label_set.push_back(make_pair(idx.spt_v[i],idx.spt_d[i]));
            }
        }
    }

    vector<vector<raw_label>> get_raw_label_list(    const vector<NodeID>&  mapper){
        raw_label_list.resize(numOfVertices);
        for(NodeID s = 0; s <numOfVertices; s++){
            const index_t_path_p &idx = index_p[s];
            _mm_prefetch(&idx.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx.spt_d[0], _MM_HINT_T0);
            _mm_prefetch(&idx.spt_p[0], _MM_HINT_T0);
            for(int i = 0; i < numOfVertices; ++i){
                //if(i>0 and idx.spt_v[i]==NULL){
                //boundary case to stop
                if(idx.spt_v[i] == numOfVertices){
                    break;
                }else{
                    raw_label_list[s].push_back(raw_label{idx.spt_v[i],idx.spt_d[i],idx.spt_p[i],-1 });
                }
            }
        }
        for(unsigned i = 0 ; i < raw_label_list.size();i ++){
            for(raw_label& r : raw_label_list[i]){
                int pre_predecessor = i;
                int predecessor = r.label_predecessor;
                while(predecessor !=  mapper[r.hub_id]){
                    const index_t_path_p &idx_from_s = index_p[predecessor];

                    _mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

                    //	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
                    for (int j = 0; ; ++j) {
                        if (idx_from_s.spt_v[j] == numOfVertices){
                            std::cout<<"Error, first node not found"<<std::endl;
                            break;
                        }
                        if (idx_from_s.spt_v[j] == r.hub_id) {
                            pre_predecessor = predecessor;
                            predecessor = idx_from_s.spt_p[j];
                            break;
                        }
                    }
                }
                r.label_first_node = pre_predecessor;
            }
        }
        return raw_label_list;
    }


    void get_label(NodeID s, vector<raw_label>& output){
        const index_t_path_p &idx = index_p[s];
        _mm_prefetch(&idx.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx.spt_p[0], _MM_HINT_T0);
        for(int i = 0; i < numOfVertices; ++i){
            //if(i>0 and idx.spt_v[i]==NULL){
            //boundary case to stop
            if(idx.spt_v[i] == numOfVertices){
                break;
            }else{
                output.push_back(raw_label{idx.spt_v[i],idx.spt_d[i],idx.spt_p[i] });
            }
        }
    }


    EdgeWeight query(NodeID s, NodeID t) {
        EdgeWeight distance = INF_WEIGHT;
        const index_t_path_p &idx_s = index_p[s];
        const index_t_path_p &idx_t = index_p[t];

        _mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

        for (int i = 0, j = 0; ; ) {
            NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

            if (v1 == numOfVertices) break;  // Sentinel

            if (v1 == v2) {
                EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
                if (td < distance) {
                    distance = td;
                    //cout << "Node: " << j << " Distance: " << idx_s.spt_d[i] << endl;
                }
                ++i;
                ++j;
            }
            else {
                i += v1 < v2 ? 1 : 0;
                j += v1 > v2 ? 1 : 0;
            }
        }
        //
        if(distance == INF_WEIGHT){
            distance = -1;
        }
        return distance;
    }



    EdgeWeight query(NodeID s, NodeID t, NodeID& meet, EdgeWeight& dis1, EdgeWeight& dis2) {
        EdgeWeight distance = INF_WEIGHT;
        vector<NodeID>& index_s = index_[s].spt_v;
        vector<EdgeWeight>& index_s_d = index_[s].spt_d;

        vector<NodeID>& index_t = index_[t].spt_v;
        vector<EdgeWeight>& index_t_d = index_[t].spt_d;
        meet = numeric_limits<NodeID>::max();
        dis1 = numeric_limits<EdgeWeight>::max();
        dis2 = numeric_limits<EdgeWeight>::max();
        for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
            if (index_s[i] == index_t[j]) {
                if (distance >(EdgeWeight)(index_s_d[i] + index_t_d[j])) {
                    distance = (EdgeWeight)(index_s_d[i] + index_t_d[j]);
                    meet = index_s[i];
                    dis1 = index_s_d[i];
                    dis2 = index_t_d[j];
                }
                ++i; ++j;
            }
            else {
                if (index_s[i] < index_t[j])
                    ++i;
                else
                    ++j;
            }
        }
        return distance;
    }


    EdgeWeight query_path(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv) {
        EdgeWeight distance = INF_WEIGHT;
        NodeID meetnode = numOfVertices;
        NodeID s_parent;
        NodeID t_parent;

        const index_t_path_p &idx_s = index_p[s];
        const index_t_path_p &idx_t = index_p[t];

        _mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

        for (int i = 0, j = 0; ; ) {
            NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

            if (v1 == numOfVertices) break;  // Sentinel

            if (v1 == v2) {
                EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
                if (td < distance) {
                    distance = td;
                    //	if (v1 < meetnode) {
                    meetnode = v1;
                    s_parent = idx_s.spt_p[i];
                    t_parent = idx_t.spt_p[j];
                    //}
                }
                ++i;
                ++j;
            }
            else {
                i += v1 < v2 ? 1 : 0;
                j += v1 > v2 ? 1 : 0;
            }
        }


        //Next, retrieve path from s - meetnode and meetnode - t.
        vector<NodeID> path_from_s;
        vector<NodeID> path_to_t;
        path_from_s.push_back(s_parent);
        path_to_t.push_back(t_parent);

        int operation = 0;

        /*	if (s == 194569 && t == 20072)
        cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

        NodeID inv_meetnode = inv[meetnode];

        while (path_from_s.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "s meet:" << path_from_s.back() << endl;*/
            const index_t_path_p &idx_from_s = index_p[path_from_s.back()];

            _mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

            //	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
            for (int i = 0; ; ++i) {
                operation++;
                if (idx_from_s.spt_v[i] == numOfVertices) break;
                if (idx_from_s.spt_v[i] == meetnode) {
                    path_from_s.push_back(idx_from_s.spt_p[i]);
                    break;
                }
            }
        }

        while (path_to_t.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "t meet:" << path_to_t.back() << endl;*/
            //	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
            const index_t_path_p &idx_to_t = index_p[path_to_t.back()];

            _mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
            for (int i = 0; ; ++i) {
                operation++;
                if (idx_to_t.spt_v[i] == numOfVertices) break;
                if (idx_to_t.spt_v[i] == meetnode) {
                    path_to_t.push_back(idx_to_t.spt_p[i]);
                    break;
                }
            }
        }

        distance = 0;
        distance += path_from_s.size() + path_to_t.size();

        return distance;

    }

    EdgeWeight query_prefix_path(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv,vector<int>&path, unsigned path_length) {


        EdgeWeight distance = INF_WEIGHT;
        NodeID meetnode = numOfVertices;
        NodeID s_parent;
        NodeID t_parent;

        const index_t_path_p &idx_s = index_p[s];
        const index_t_path_p &idx_t = index_p[t];

        _mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

        for (int i = 0, j = 0; ; ) {
            NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

            if (v1 == numOfVertices) break;  // Sentinel

            if (v1 == v2) {
                EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
                if (td < distance) {
                    distance = td;
                    //	if (v1 < meetnode) {
                    meetnode = v1;
                    s_parent = idx_s.spt_p[i];
                    t_parent = idx_t.spt_p[j];
                    //}
                }
                ++i;
                ++j;
            }
            else {
                i += v1 < v2 ? 1 : 0;
                j += v1 > v2 ? 1 : 0;
            }
        }


        //Next, retrieve path from s - meetnode and meetnode - t.


        int operation = 0;

        /*	if (s == 194569 && t == 20072)
        cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

        NodeID inv_meetnode = inv[meetnode];
//        vector<NodeID> path_from_s;
        vector<NodeID> path_to_t;
        if(s !=inv_meetnode)
            path.push_back(s);

        path.push_back(s_parent);
        if(path.size() == path_length){
            return distance;
        }

        while (path.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "s meet:" << path_from_s.back() << endl;*/
            const index_t_path_p &idx_from_s = index_p[path.back()];

            _mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

            //	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
            for (int i = 0; ; ++i) {
                operation++;
                if (idx_from_s.spt_v[i] == numOfVertices) break;
                if (idx_from_s.spt_v[i] == meetnode) {
                    path.push_back(idx_from_s.spt_p[i]);
                    if(path.size() == path_length){
                        return distance;
                    }
                    break;
                }
            }
        }






        if (t != inv_meetnode)
            path_to_t.push_back(t);

        path_to_t.push_back(t_parent);

        while (path_to_t.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "t meet:" << path_to_t.back() << endl;*/
            //	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
            const index_t_path_p &idx_to_t = index_p[path_to_t.back()];

            _mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
            for (int i = 0; ; ++i) {
                operation++;
                if (idx_to_t.spt_v[i] == numOfVertices) break;
                if (idx_to_t.spt_v[i] == meetnode) {
                    path_to_t.push_back(idx_to_t.spt_p[i]);
                    break;
                }
            }
        }



        std::reverse(path_to_t.begin(),path_to_t.end());
        unsigned tmp_size = path.size();
        if(path_length > path.size()+path_to_t.size() -1){
            path.resize( path.size()+path_to_t.size() -1);
            std::copy(path_to_t.begin()+1,path_to_t.end(),path.begin()+tmp_size);

        }else{
            path.resize( path_length);
            std::copy(path_to_t.begin()+1,path_to_t.begin()+1 + path_length - tmp_size,path.begin()+tmp_size);

        }




        return distance;
    }

    EdgeWeight query_path(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv,vector<int>&path) {


        EdgeWeight distance = INF_WEIGHT;
        NodeID meetnode = numOfVertices;
        NodeID s_parent;
        NodeID t_parent;

        const index_t_path_p &idx_s = index_p[s];
        const index_t_path_p &idx_t = index_p[t];

        _mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

        for (int i = 0, j = 0; ; ) {
            NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

            if (v1 == numOfVertices) break;  // Sentinel

            if (v1 == v2) {
                EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
                if (td < distance) {
                    distance = td;
                    //	if (v1 < meetnode) {
                    meetnode = v1;
                    s_parent = idx_s.spt_p[i];
                    t_parent = idx_t.spt_p[j];
                    //}
                }
                ++i;
                ++j;
            }
            else {
                i += v1 < v2 ? 1 : 0;
                j += v1 > v2 ? 1 : 0;
            }
        }


        //Next, retrieve path from s - meetnode and meetnode - t.

        //check whether it's reachable
        if(distance == INF_WEIGHT){
            distance = -1;
            path.push_back(-1);
            return distance;
        }
        int operation = 0;

        /*	if (s == 194569 && t == 20072)
        cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

        NodeID inv_meetnode = inv[meetnode];
//        vector<NodeID> path_from_s;
        vector<NodeID> path_to_t;
        if(s !=inv_meetnode)
            path.push_back(s);
        path.push_back(s_parent);


        while (path.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "s meet:" << path_from_s.back() << endl;*/
            const index_t_path_p &idx_from_s = index_p[path.back()];

            _mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

            //	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
            for (int i = 0; ; ++i) {
                operation++;
                if (idx_from_s.spt_v[i] == numOfVertices) break;
                if (idx_from_s.spt_v[i] == meetnode) {
                    path.push_back(idx_from_s.spt_p[i]);
                    break;
                }
            }
        }






        if (t != inv_meetnode)
            path_to_t.push_back(t);

        path_to_t.push_back(t_parent);

        while (path_to_t.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "t meet:" << path_to_t.back() << endl;*/
            //	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
            const index_t_path_p &idx_to_t = index_p[path_to_t.back()];

            _mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
            for (int i = 0; ; ++i) {
                operation++;
                if (idx_to_t.spt_v[i] == numOfVertices) break;
                if (idx_to_t.spt_v[i] == meetnode) {
                    path_to_t.push_back(idx_to_t.spt_p[i]);
                    break;
                }
            }
        }


//        distance = 0;
//        distance += path_from_s.size() + path_to_t.size();

        std::reverse(path_to_t.begin(),path_to_t.end());
        unsigned tmp_size = path.size();
        path.resize( path.size()+path_to_t.size() -1);
//        std::copy(path_from_s.begin(),path_from_s.end()-1,path.begin());
        std::copy(path_to_t.begin()+1,path_to_t.end(),path.begin()+tmp_size);
//		return distance;


        return distance;

    }






    void query_path(NodeID s,NodeID s_parent, NodeID t,NodeID t_parent, NodeID meetnode, const vector<NodeID>& rank, const vector<NodeID>& inv, const vector<polyanya::Point>& turningpoint ,vector<polyanya::Point>&  shortestPath) {

        NodeID inv_meetnode = inv[meetnode];
//        vector<NodeID> path_from_s;
        vector<NodeID> path_to_t;
        if(s !=inv_meetnode){
            shortestPath.push_back(turningpoint[s]);
            shortestPath.back().vertexId = s;
        }
        shortestPath.push_back(turningpoint[s_parent]);
        shortestPath.back().vertexId = s_parent;

        while (shortestPath.back().vertexId != inv_meetnode) {
            const index_t_path_p &idx_from_s = index_p[shortestPath.back().vertexId];

            _mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

            //	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
            for (int i = 0; ; ++i) {
                if (idx_from_s.spt_v[i] == numOfVertices) break;
                if (idx_from_s.spt_v[i] == meetnode) {
                    shortestPath.push_back(turningpoint[idx_from_s.spt_p[i]]);
                    shortestPath.back().vertexId = idx_from_s.spt_p[i];
                    break;
                }
            }
        }

        if (t != inv_meetnode)
            path_to_t.push_back(t);

        path_to_t.push_back(t_parent);

        while (path_to_t.back() != inv_meetnode) {
            const index_t_path_p &idx_to_t = index_p[path_to_t.back()];
            _mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
            for (int i = 0; ; ++i) {
                if (idx_to_t.spt_v[i] == numOfVertices) break;
                if (idx_to_t.spt_v[i] == meetnode) {
                    path_to_t.push_back(idx_to_t.spt_p[i]);
                    break;
                }
            }
        }

        std::reverse(path_to_t.begin(),path_to_t.end());
        for(unsigned int i = 1; i < path_to_t.size(); i ++){
            shortestPath.push_back(turningpoint[path_to_t[i]]);
            shortestPath.back().vertexId =path_to_t[i];
        }

    }


    EdgeWeight query_path_check(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv) {


        EdgeWeight distance = INF_WEIGHT;
        NodeID meetnode = numOfVertices;
        NodeID s_parent;
        NodeID t_parent;

        const index_t_path_p &idx_s = index_p[s];
        const index_t_path_p &idx_t = index_p[t];

        _mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
        _mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

        for (int i = 0, j = 0; ; ) {
            NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

            if (v1 == numOfVertices) break;  // Sentinel

            if (v1 == v2) {
                EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
                if (td < distance) {
                    distance = td;
                    //	if (v1 < meetnode) {
                    meetnode = v1;
                    s_parent = idx_s.spt_p[i];
                    t_parent = idx_t.spt_p[j];
                    //}
                }
                ++i;
                ++j;
            }
            else {
                i += v1 < v2 ? 1 : 0;
                j += v1 > v2 ? 1 : 0;
            }
        }

        NodeID inv_meetnode = inv[meetnode];
        //Next, retrieve path from s - meetnode and meetnode - t.
        vector<NodeID> path_from_s;
        vector<NodeID> path_to_t;
        if(s !=inv_meetnode)
            path_from_s.push_back(s);
        path_from_s.push_back(s_parent);
        if (t != inv_meetnode)
            path_to_t.push_back(t);

        path_to_t.push_back(t_parent);

        /*	if (s == 194569 && t == 20072)
        cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/


        while (path_from_s.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "s meet:" << path_from_s.back() << endl;*/
            const index_t_path_p &idx_from_s = index_p[path_from_s.back()];

            _mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

            //	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
            for (int i = 0; ; ++i) {
                if (idx_from_s.spt_v[i] == numOfVertices) break;
                if (idx_from_s.spt_v[i] == meetnode) {
                    path_from_s.push_back(idx_from_s.spt_p[i]);
                    break;
                }
            }
        }

        while (path_to_t.back() != inv_meetnode) {
            /*if (s == 194569 && t == 20072)
            cout << "t meet:" << path_to_t.back() << endl;*/
            //	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
            const index_t_path_p &idx_to_t = index_p[path_to_t.back()];
            _mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
            _mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
            for (int i = 0; ; ++i) {
                if (idx_to_t.spt_v[i] == numOfVertices) break;
                if (idx_to_t.spt_v[i] == meetnode) {
                    path_to_t.push_back(idx_to_t.spt_p[i]);
                    break;
                }
            }
        }

        //return distance;

        EdgeWeight alldis = 0;

        if (path_from_s.size() == 1)
            if (s != inv_meetnode)
                alldis += query_p(s, inv_meetnode);

        if (path_to_t.size() == 1)
            if (t != inv_meetnode)
                alldis += query_p(t, inv_meetnode);

        for (int i = 0; i < path_from_s.size() - 1; ++i) {
            alldis += query_p(path_from_s[i], path_from_s[i + 1]);
            //cout << "s " << path_from_s[i] << "," << path_from_s[i + 1] << endl;
        }
        for (int i = 0; i < path_to_t.size() - 1; ++i) {
            alldis += query_p(path_to_t[i], path_to_t[i + 1]);

            //cout <<"t " <<  path_to_t[i] << "," << path_to_t[i + 1] << endl;
        }
        /*if (distance != alldis)
            cout << "a?" << endl;*/
        //cout << distance << "," << alldis << "," << path_from_s.size() + path_to_t.size() << endl;
//		cout << s << "," << t << "," << inv_meetnode << "   " << distance << "vs." << alldis << endl;

        return distance;
    }


    double avg_size() {
        double total = 0;
        for (int i = 0; i < numOfVertices; ++i) total += index_[i].spt_v.size();

        double avg = total / numOfVertices - 1; // We do not count the trivial label (V, INF_WEIGHT).

        return avg;
    }
    double total_size() {
        double total = 0;
        for (int i = 0; i < numOfVertices; ++i) total += index_[i].spt_v.size();
        return total;
    }
    /*
    NodeID max_size() {
    NodeID maxsize = numeric_limits<NodeID>::min();
    for (int i = 0; i < V; ++i) maxsize = max(maxsize, index_[i].spt_v.size());
    return maxsize;
    }*/

    void append(NodeID v, NodeID root, EdgeWeight distance) {
        index_[v].spt_v.push_back(root);
        index_[v].spt_d.push_back(distance);
    }

    void print_stat() {
        cout << "Average Label Size: " << avg_size() << endl;
        //cout << "Maximum Label Size: " << max_size() << endl;
    }

    void Free() {
        if (index_.size() == 0) return;
        for (int v = 0; v < numOfVertices; ++v) {
            index_[v].spt_v.clear();
            index_[v].spt_d.clear();
        }
        index_.clear();
    }

    void save_labels(const char* save_filename) {
        ofstream ofs(save_filename, ios::binary | ios::out);

        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        for (NodeID v = 0; v < numOfVertices; ++v) {
            NodeID isize = index_[v].size();
            ofs.write((const char*)&isize, sizeof(isize));
            for (NodeID i = 0; i < index_[v].size(); ++i) {
                ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
                ofs.write((const char*)&index_[v].spt_p[i], sizeof(index_[v].spt_p[i]));
                ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
            }
        }
        ofs.close();
    }
    //label memory with calculate euclidean symbol
    void save_labels(const char* save_filename,  double& average_label, unsigned long long& memory, double& max_label, vector<vector<int>> &visible_vertex, double &repeated_label) {
        ofstream ofs(save_filename, ios::binary | ios::out);
        unsigned long long label_size = 0;
        unsigned long long temp_label_size = 0;
        unsigned long long covisible_label = 0;
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        for (NodeID v = 0; v < numOfVertices; ++v) {
            NodeID isize = index_[v].size();
            if(isize > temp_label_size){
                temp_label_size = isize;
            }
            ofs.write((const char*)&isize, sizeof(isize));
            //memory doesn't include predecessor as it's distance retrieval
            for (NodeID i = 0; i < index_[v].size(); ++i) {
                ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
                ofs.write((const char*)&index_[v].spt_p[i], sizeof(index_[v].spt_p[i]));
                ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
                label_size ++;
                //cout << "Node ID: " << v << " hub: " << index_[v].spt_v[i] << " predecessor: " << index_[v].spt_p[i]<<
                //" distance: " << index_[v].spt_d[i] << endl;
                //calculate euclidean symbol using visibility graph
                for (int j = 0; j < visible_vertex[v].size(); ++j){
                    if(index_[v].spt_v[i] == visible_vertex[v][j]){
                        covisible_label ++;
                    }
                }
                memory= memory + sizeof(index_[v].spt_v[i])  + sizeof(index_[v].spt_d[i]);

            }
        }
        repeated_label = covisible_label/(double)label_size;
        cout << "Total Label: " << label_size << endl;
        cout << "Covisible Label: " << covisible_label << endl;
        average_label = label_size/(double)numOfVertices;
        max_label = temp_label_size;
        ofs.close();
    }
    //path memory
    void save_labels(const char* save_filename,  double& average_label, unsigned long long& path_memory, unsigned long long& original_memory) {
        ofstream ofs(save_filename, ios::binary | ios::out);
        unsigned long long label_size = 0;
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        for (NodeID v = 0; v < numOfVertices; ++v) {
            NodeID isize = index_[v].size();
            ofs.write((const char*)&isize, sizeof(isize));
            for (NodeID i = 0; i < index_[v].size(); ++i) {
                ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
                ofs.write((const char*)&index_[v].spt_p[i], sizeof(index_[v].spt_p[i]));
                ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
                label_size ++;
                if(index_[v].spt_v[i] != numOfVertices) {
                    path_memory = path_memory + sizeof(index_[v].spt_v[i]) + sizeof(index_[v].spt_p[i]) + sizeof(index_[v].spt_d[i]);
                    original_memory = original_memory + sizeof(index_[v].spt_v[i]) + sizeof(index_[v].spt_d[i]);
                }
            }
        }
        average_label = avg_size();
        std::cout<<"Path memory: "<< path_memory<<" Bytes"<<std::endl;
        std::cout<<"Original memory: "<< original_memory<<" Bytes"<<std::endl;
        ofs.close();
    }

    void load_labels(const char* load_filename) {
        /*	for (NodeID v = 0; v < numOfVertices; ++v) {
        free(index_p[v].spt_v);
        free(index_p[v].spt_d);
        }
        */
        //free(index_p);
        index_p = NULL;

        ifstream ifs(load_filename);
        if(!ifs){
            perror("Failed because: ");
            exit(1);
        }
        NodeID isize = 0;
        ifs.read((char*)&isize, sizeof(isize));
        numOfVertices = isize;
        //numOfVertices = 65;
//        std::cout<<numOfVertices<<std::endl;
        index_p = (index_t_path_p*)memalign(64, numOfVertices * sizeof(index_t_path_p));

        for (NodeID v = 0; v < numOfVertices; ++v) {
            index_t_path_p &idx = index_p[v];
            ifs.read((char*)&isize, sizeof(isize));

            idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
            idx.spt_p = (NodeID*)memalign(64, isize * sizeof(NodeID));
            idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

            //	index_[v].spt_v.resize(isize);
            //	index_[v].spt_d.resize(isize);

            for (NodeID i = 0; i < isize; ++i) {
                NodeID hub;
                NodeID hub_parent;
                EdgeWeight hub_weight;
                ifs.read((char*)&hub, sizeof(hub));
                ifs.read((char*)&hub_parent, sizeof(hub_parent));
                ifs.read((char*)&hub_weight, sizeof(hub_weight));
                //index_[v].spt_v[i] = hub;
                //index_[v].spt_d[i] = hub_weight;
                idx.spt_v[i] = hub;
                idx.spt_p[i] = hub_parent;
                idx.spt_d[i] = hub_weight;
            }
        }
        ifs.close();


    }

    void save_labels_iteration_stats(const char* save_filename) {

        vector<NodeID> stat(numOfVertices);
        for (NodeID v = 0; v < numOfVertices; ++v) {
            for (NodeID i = 0; i < index_[v].size(); ++i)
                stat[index_[v].spt_v[i]]++;
        }

        ofstream ofs(save_filename);

        for (NodeID v = 0; v < numOfVertices; ++v) {
            ofs << stat[v] << endl;
        }
        ofs.close();
    }

    EdgeWeight query_with_info(NodeID s, NodeID t, query_info& q_info) {

        double stime = GetCurrentTimeSec();

        EdgeWeight distance = INF_WEIGHT;
        vector<NodeID>& index_s = index_[s].spt_v;
        vector<EdgeWeight>& index_s_d = index_[s].spt_d;

        vector<NodeID>& index_t = index_[t].spt_v;
        vector<EdgeWeight>& index_t_d = index_[t].spt_d;

        q_info.meet_node = numOfVertices;
        double meet_distance;

        for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
            if (index_s[i] == index_t[j]) {
                meet_distance = (EdgeWeight)(index_s_d[i++] + index_t_d[j++]);
                if (distance >  meet_distance) {
                    distance = meet_distance;
                    q_info.meet_node = index_s[i];
                }
            }
            else {
                if (index_s[i] < index_t[j])
                    ++i;
                else
                    ++j;
            }
        };

        stime = GetCurrentTimeSec() - stime;

        q_info.time_cost = stime;

        if (index_s.size() < index_t.size())
            q_info.search_len = index_s.size();
        else
            q_info.search_len = index_t.size();

        return distance;
    }

    void label_memory(unsigned long long& path_memory) {
        for (NodeID v = 0; v < numOfVertices; ++v) {
            NodeID isize = index_[v].size();
            for (NodeID i = 0; i < index_[v].size(); ++i) {
                if(index_[v].spt_v[i] != numOfVertices) {
                    path_memory = path_memory + sizeof(index_[v].spt_v[i]) + sizeof(index_[v].spt_p[i]) + sizeof(index_[v].spt_d[i]);
                }
            }
        }
    }


};

#endif