#!/bin/bash
echo Creating directory for visibility graph, hub label and EHL...
mkdir dataset/visibility_graph
mkdir dataset/hub_label
mkdir dataset/ehl
mkdir dataset/visibility_graph/$1
mkdir dataset/hub_label/$1
mkdir dataset/ehl/$1
mkdir dataset/ehl/$1/workload

for file in dataset/merged-mesh/$1/* ;do
  echo Creating visibility graph and hub labels for $file ..........
  array=(${file//// });
  directory_name=${array[2]};
  f="$(basename -- $file)";
  array2=(${f//-/ });
  map_name=${array2[0]};

./bin/build_visibility_graph dataset/merged-mesh/$directory_name/$map_name-merged.mesh dataset/grid/$directory_name/$map_name.map dataset/visibility_graph/$directory_name/$map_name.vis
./bin/hub_label_converter dataset/visibility_graph/$directory_name/$map_name.vis dataset/visibility_graph/$directory_name/$map_name.txt
./bin/construct_hl dataset/visibility_graph/$directory_name/$map_name.txt dataset/hub_label/$directory_name/$map_name dataset/hub_label/result/$directory_name/$map_name $map_name $directory_name


  ./bin/build_grid_based_hub_labelling $directory_name $map_name 1;
  #Input: benchmark map grid_size(to determine preprocessed EHL) memory_target% cluster-x
  ./bin/build_query_workload_index $directory_name $map_name 1 80 2


  echo Finish preprocessing .......................................................................

done