#!/bin/bash

for file in dataset/merged-mesh/$1/* ;do
  echo Processing query for $file ..........
  array=(${file//// });
  directory_name=${array[2]};
  f="$(basename -- $file)";
  array2=(${f//-/ });
  map_name=${array2[0]};
  echo $map_name

 #Input: benchmark map grid_size(to determine preprocessed EHL) cluster-x memory_target%
./bin/testEHL $directory_name $map_name 1 2 80

done