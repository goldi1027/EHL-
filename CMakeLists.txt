project(EHL)
cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-pedantic -Wall -Wno-strict-aliasing -Wno-long-long -Wno-deprecated -Wno-deprecated-declarations")
#set(DEV_CXX_FLAGS "-g -ggdb -O0 -fno-omit-frame-pointer -fsanitize=address")
set(DEV_CXX_FLAGS "-g -ggdb -O0 -fno-omit-frame-pointer")
set(RELEASE_CXX_FLAGS "-O3 -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${DEV_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${RELEASE_CXX_FLAGS}")


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/../bin)


if(UNIX AND NOT APPLE)
    message(STATUS ">>> Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include -fopenmp")
    find_library(omp_lib omp "/usr/local/opt")
    include_directories("./boost_1_64_0")

endif ()


if(APPLE)
    include_directories("./boost_1_64_0")
    include_directories("/usr/local/Cellar/google-sparsehash/2.0.4/include")
    include_directories("/usr/local/Cellar/libomp/14.0.6/include")
    find_library(omp_lib omp "/usr/local/opt")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
endif()


include_directories(helpers)
include_directories(search)
include_directories(structs)
include_directories(hub_label_structs)

set(SRC
        helpers/cfg.cpp
        helpers/cfg.h
        helpers/cpool.h
        helpers/geometry.cpp
        helpers/geometry.h
        helpers/rtree.h
        helpers/scenario.cpp
        helpers/scenario.h
        helpers/timer.cpp
        helpers/timer.h
        helpers/timestamp_flag.h
        helpers/constants.h
        helpers/id_queue.h
        helpers/helper.h

        search/expansion.cpp
        search/expansion.h
        search/searchinstance.cpp
        search/searchinstance.h
        search/visibleSearchInstance.cpp
        search/visibleSearchInstance.h
        search/dijkstra.h
        search/visibleAreaSearchInstance.h
        search/visibleAreaSearchInstance.cpp
        search/ebhl_query_v2.h
        search/query_workload_query.h

        structs/vec_io.h
        structs/consts.h
        structs/mesh.cpp
        structs/mesh.h
        structs/point.h
        structs/polygon.h
        structs/searchnode.h
        structs/successor.h
        structs/vertex.h
        structs/graph.h
        structs/graph.cpp
        structs/binary_search.h
        structs/range.h
        structs/grid_label.h
        structs/ebhl.h
        structs/edge.h
        structs/query_workload.h

        hub_label_structs/coverage_ordering_path.h
        hub_label_structs/graph.h
        hub_label_structs/graph_search.h
        hub_label_structs/heap.h
        hub_label_structs/labels.h
        hub_label_structs/paras.h
        hub_label_structs/time_util.h

        )

add_executable(build_visibility_graph ${SRC} build_visibility_graph.cpp)
add_executable(construct_hl ${SRC} construct_hl.cpp)
add_executable(hub_label_converter ${SRC} hub_label_converter.cpp)
add_executable(testEHL ${SRC} testEHL.cpp)
add_executable(build_grid_based_hub_labelling ${SRC} build_grid_based_hub_labelling.cpp)
add_executable(build_query_workload_index ${SRC} build_query_workload_index.cpp)

if (omp_lib)
    target_link_libraries(build_grid_based_hub_labelling ${omp_lib})
    target_link_libraries(build_visibility_graph ${omp_lib})
    target_link_libraries(testEHL ${omp_lib})
endif()
