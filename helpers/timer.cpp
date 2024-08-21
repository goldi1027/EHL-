/*
 Compromise-free Pathfinding on a Navigation Mesh
 Authors: Michael Cui, Daniel Harabor and Alban Grastien
 Published venue: Proceedings of the Twenty-Sixth International Joint Conference on Artificial Intelligence, 2017
 Link to source code: https://bitbucket.org/dharabor/pathfinding/src/master/anyangle/polyanya/

 This implementation of Polyanya is licensed under MIT.
 Several source files from Daniel Harabor's Warthog project were used this project - these files are also licensed under MIT. These files are: helpers/cfg.cpp, helpers/cfg.h, helpers/cpool.h, helpers/timer.cpp and helpers/timer.h.
 */

#include <stdio.h>
//#include "timer.h"
#include <chrono>
#include "timer.h"

warthog::timer::timer()
{

//#ifdef OS_MAC
//    start_time = stop_time = 0;
//    mach_timebase_info(&timebase);
//
//#else
////    start_time.tv_sec = 0;
////    start_time.tv_nsec = 0;
////    stop_time.tv_sec = 0;
////    stop_time.tv_nsec = 0;
//#endif

}

double
warthog::timer::get_time_nano()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();
//#ifdef OS_MAC
//    uint64_t raw_time = mach_absolute_time();
//    return (double)(raw_time * timebase.numer / timebase.denom);
//#else
//    timespec raw_time;
//    clock_gettime(CLOCK_MONOTONIC , &raw_time);
//    return (double)(raw_time.tv_nsec);
//#endif
}

void warthog::timer::start()
{
    start_time = std::chrono::steady_clock::now();
    stop_time = start_time;
//#ifdef OS_MAC
//    start_time = mach_absolute_time();
//    stop_time = start_time;
//#else
////    clock_gettime(CLOCK_MONOTONIC , &start_time);
//#endif
}

void warthog::timer::stop()
{
    stop_time =  std::chrono::steady_clock::now();
//#ifdef OS_MAC
//    stop_time = mach_absolute_time();
//#else
////    clock_gettime(CLOCK_MONOTONIC , &stop_time);
////    auto stime = std::chrono::steady_clock::now();
////    auto etime = std::chrono::steady_clock::now();
//#endif
}

double warthog::timer::current_time_nano(std::chrono::time_point<std::chrono::steady_clock> current_time){
    return std::chrono::duration_cast<std::chrono::nanoseconds>(current_time - start_time).count();
}
double warthog::timer::elapsed_time_nano()
{

    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();

//#ifdef OS_MAC
//    uint64_t elapsed_time = stop_time - start_time;
//    return (double)(elapsed_time * timebase.numer / timebase.denom);
////    Nanoseconds nanosecs = AbsoluteToNanoseconds(*(AbsoluteTime*)&elapsed_time);
////    return (double) UnsignedWideToUInt64(nanosecs) ;
//
//#else
////    if ((stop_time.tv_nsec-start_time.tv_nsec)<0)
////        return (double)(1000000000+stop_time.tv_nsec-start_time.tv_nsec);
////    else
////        return (double)(stop_time.tv_nsec - start_time.tv_nsec);
//#endif
}

void warthog::timer::reset()
{
//#ifdef OS_MAC
//    start_time = stop_time = 0;
//#else
////    start_time.tv_sec = 0;
////    start_time.tv_nsec = 0;
////    stop_time.tv_sec = 0;
////    stop_time.tv_nsec = 0;
//#endif
}

double
warthog::timer::elapsed_time_micro()
{
    return elapsed_time_nano() / 1000.0;
}

