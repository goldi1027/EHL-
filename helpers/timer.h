// timer.h
//
// A cross-platform monotonic wallclock timer.
// Currently supports nanoseconds resolution.
//
// Reference doco for timers on OSX:
// https://developer.apple.com/library/mac/qa/qa1398/_index.html
// https://developer.apple.com/library/mac/technotes/tn2169/_index.html#//apple_ref/doc/uid/DTS40013172-CH1-TNTAG5000
//

#pragma once
#include <chrono>
//#ifdef OS_MAC
////#include <CoreServices/CoreServices.h>
//#include <mach/mach.h>
//#include <mach/mach_time.h>
//
//#else
////#include <time.h>
//#endif

namespace warthog
{

class timer
{
   std::chrono::time_point<std::chrono::steady_clock> stop_time;
   std::chrono::time_point<std::chrono::steady_clock> start_time;
//#ifdef OS_MAC
//    uint64_t start_time;
//    uint64_t stop_time;
//    mach_timebase_info_data_t timebase;
//#else
//    static std::chrono::time_point<std::chrono::steady_clock> stop_time;
//    static std::chrono::time_point<std::chrono::steady_clock> start_time;
////    timespec stop_time;
////    timespec start_time;
//#endif

public:
    timer();
    void reset();
    void start();
    void stop();

    double current_time_nano(std::chrono::time_point<std::chrono::steady_clock> current_time);
    double elapsed_time_nano();
    double elapsed_time_micro();
    double get_time_nano();
};

}
