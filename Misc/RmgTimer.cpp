
#include "RmgTimer.h"
#include <string>

#if 1
volatile double start_time;
volatile double end_time;
const char *sname;
std::unordered_map<std::string, double> RmgTimer::timings[MAX_RMG_THREADS+1];

RmgTimer::RmgTimer(const char *fname) {
    struct timeval t1;
    gettimeofday (&t1, NULL);
    sname = fname;
    start_time = t1.tv_sec + 1e-6 * t1.tv_usec;
}

RmgTimer::~RmgTimer(void) {
    int tid;
    BaseThread T(0), *Tptr;
    struct timeval t1;
    gettimeofday (&t1, NULL);
    end_time = t1.tv_sec + 1e-6 * t1.tv_usec;
    Tptr = T.get_thread_control();

    // if Tptr is null that means we are being called from the main program
    if(!Tptr) {
        tid = MAX_RMG_THREADS;
    }
    else {
        tid = T.get_thread_tid();
        if(tid < 0) tid=MAX_RMG_THREADS;
    }

    if(timings[tid].count(sname)) {
        RmgTimer::timings[tid][sname] += end_time - start_time;
    }
    else {
        RmgTimer::timings[tid][sname] = end_time - start_time;
    }
}

// Iterates over list and prints per thread timing statistics
void RmgTimer::PrintTimings(void) {

    for(auto it = RmgTimer::timings[MAX_RMG_THREADS].cbegin(); it != RmgTimer::timings[MAX_RMG_THREADS].cend(); ++it)
    {
        cout << "  Main thread time for " << it->first << " = " << it->second << "\n";
    }

}


extern "C" void CompatRmgTimerPrint(void) {
    RmgTimer RT("Print timings");
    RT.PrintTimings();
}

#endif

