
#include <sys/time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <functional>
#include <string>
#include <iomanip>
#include "RmgTimer.h"

using namespace std;

volatile double start_time;
volatile double end_time;
const char *sname;
std::unordered_map<std::string, double> RmgTimer::timings[MAX_RMG_THREADS+1];



RmgTimer::RmgTimer(const char *what) 
{
    struct timeval t1;
    gettimeofday (&t1, NULL);
    sname = what;
    start_time = t1.tv_sec + 1e-6 * t1.tv_usec;
}

RmgTimer::~RmgTimer(void) 
{
    if(!sname) return;  // Bit of a hack for the temporary C interfaces
    struct timeval t1;
    gettimeofday (&t1, NULL);
    end_time = t1.tv_sec + 1e-6 * t1.tv_usec;

    int tid;
    BaseThread T(0), *Tptr;
    Tptr = T.get_thread_control();

    // if Tptr is null that means we are being called from the main program
    if(!Tptr) {
        tid = 0;
    }
    else {
        tid = T.get_thread_tid();
        if(tid < 0) {
            tid = 0;
        }
        else {
            tid++;      // Increment by 1 since we store the main thread values in the zeroth slot
        }
    }

    if(timings[tid].count(sname)) {
        RmgTimer::timings[tid][sname] += end_time - start_time;
    }
    else {
        RmgTimer::timings[tid][sname] = end_time - start_time;
    }
}

unordered_map<std::string, double> *RmgTimer::get_map(void)
{
    return RmgTimer::timings;
}

// Temporary until C to C++ migration is completed. Use carefully!
extern "C" void *BeginRmgTimer(const char *what) 
{
    RmgTimer *RT = new RmgTimer(what);
    return (void *)RT;
}
extern "C" void EndRmgTimer(void *ptr) 
{
    RmgTimer *RT = (RmgTimer *)ptr;
    delete(RT);
}

