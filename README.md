# Euclidean Hub Labeling Star (EHL*)

## Introduction

Implementation of Euclidean Hub Labeling Star (EHL*). EHL is an ultrafast optimal path finding algorithm used for finding shortest paths in Euclidean space.
EHL utilises visibility graph and extends hub labeling[1] to build the index used for online pathfinding. The algorithm and details of EHL is available in [2]. 
EHL* is an improved version of EHL that offers the flexibility to build the index within a specified memory budget while optimizing query runtime.
EHL code is accessed and retrieved from https://github.com/goldi1027/EHL. 


## Dataset

The four benchmarks used by EHL (dao, da2, sc1, and bgmaps) are retrieved from MovingAI (https://movingai.com).

The synthetic query sets are generated based on the four benchmarks mentioned above. Due to size limit of supplementary material, we have uploaded the full synthetic dataset of our cluster data at https://drive.google.com/file/d/1BZ3f9BaY4O1hf8n2dR4PI3INU5plvgpt/view?usp=sharing

The merged-meshes are provided by the authors of Polyanya [3] and available from repository (https://bitbucket.org/%7B3c286763-d509-45c2-b036-75814ce8955f%7D/)
The scenarios used for testing and experiments in the dataset to reproduce our results are also provided in the code.
## Requirements
The following libraries need to be installed in order to reproduce the implementation results.
For installation of the libraries, please follow the guidelines provided by the links.

- CMake: https://cmake.org
- OpenMP: https://www.openmp.org
- Google Sparsehash: https://github.com/justinsb/google-sparsehash
- Boost geometry - EHL uses boost geometry library to find and check visibility area as required in the algorithm. However, the latest version doesn't handle some precise cases and only version 1.64 is accurate to our precision needs. Boost version 1.64 can be downloaded from here (https://www.boost.org/users/history/version_1_64_0.html)

After installing the libraries, as EHL* uses Cmake to compile, certain changes and modifications should be made to CMakeLists.txt based on different machine settings.

Then with the Makefile which we have provided, the code could be compiled using "make fast".

## Running

Currently, we provide three bash scripts to quickly reproduce the experimental results reported in our paper.

bash preprocessing_EHL*.sh [BENCHMARK_SUITE] 
e.g., run "bash preprocessing_EHL*.sh dao" This bash command creates all the indexes (visibility graph, hub label, EHL and EHL*) needed for EHL* for all the maps in the benchmark suite (dao).

bash benchmark_EHL*.sh [BENCHMARK_SUITE] 
e.g., run "bash benchmark_EHL.sh dao" This bash command runs queries for EHL* for all the maps in the benchmark suite (dao) using the queries from synthetic cluster data. The output shows the average runtime of an entire map's queries to the console. 

bash clean_index.sh [BENCHMARK_SUITE]
e.g., run "bash clean_index.sh". This bash command deletes and cleans all indexes of EHL* directories.


## References
[1] Ye Li, Leong Hou U, Man Lung Yiu, Ngai Meng Kou: An Experimental Study on Hub Labeling based Shortest Path Algorithms. Proc. VLDB Endow. 11(4): 445-457 (2017)

[2] Jinchun Du, Bojie Shen, Muhammad Aamir Cheema: Ultrafast Euclidean Shortest Path Computation using Hub Labeling. AAAI 2023

[3] Michael Cui, Daniel Damir Harabor, Alban Grastien: Compromise-free Pathfinding on a Navigation Mesh. IJCAI 2017: 496-502
