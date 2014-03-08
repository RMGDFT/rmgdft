
#include "RmgTimer.h"
#include "BaseGrid.h"
#include <string>

volatile double start_time;
volatile double end_time;
const char *sname;
std::unordered_map<std::string, double> RmgTimer::timings[MAX_RMG_THREADS+1];



RmgTimer::RmgTimer(const char *fname) 
{
    struct timeval t1;
    gettimeofday (&t1, NULL);
    sname = fname;
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


// Iterates over list and prints per thread timing statistics
void RmgTimer::PrintTimings(void) {

    int tid;
    BaseThread T(0);
    BaseGrid G;

    if(G.get_gridpe() == 0) {

        for(tid=0;tid <= T.get_threads_per_node();tid++) {

            for(auto it = RmgTimer::timings[tid].cbegin(); it != RmgTimer::timings[tid].cend(); ++it)
            {
                cout << "  Time spent by thread " << tid << " in " << it->first << " = " << it->second << "\n";
            }

        }

    }
}

extern "C" void CompatRmgTimerPrint(void) {
    RmgTimer RT("Print timings");
    RT.PrintTimings();
}

