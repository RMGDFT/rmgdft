
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


// Iterates over list and prints per thread timing statistics
void RmgTimer::PrintTimings(const char *outfile) {

    int tid;
    BaseThread T(0);
    BaseGrid G;
    std::ofstream logfile;
    logfile.open(outfile, std::ofstream::out | std::ofstream::app);

    if(G.get_gridpe() == 0) {

        // Have to do some manipulations to compute thread min/max/average and get things properly sorted
        std::map <std::string, double> tmain;
        std::map <std::string, double> tmin;
        std::map <std::string, double> tmax;
        std::map <std::string, double> tavg;

        // Reset main thread (tid=0) into an ordered map
        for(auto it = RmgTimer::timings[0].cbegin(); it != RmgTimer::timings[0].cend(); ++it) {
            tmain[it->first] = it->second;
        }

        size_t maxlen = 0;
        for(auto it = RmgTimer::timings[1].cbegin(); it != RmgTimer::timings[1].cend(); ++it) {
            tmin[it->first] = 99999999999.0;
            tmax[it->first] = 0.0;
            tavg[it->first] = 0.0;
            if(maxlen < it->first.length()) maxlen = it->first.length();
        }
        maxlen += 2; 

        for(tid=1;tid <= T.get_threads_per_node();tid++) {
            for(auto it = RmgTimer::timings[tid].cbegin(); it != RmgTimer::timings[tid].cend(); ++it) {
                if(it->second < tmin[it->first]) tmin[it->first] = it->second;
                if(it->second > tmax[it->first]) tmax[it->first] = it->second;
                tavg[it->first] += it->second;
            }
        }

        logfile << "\n\n";
        logfile << std::fixed << std::setprecision(2);
        logfile << "------------------------- TIMING INFORMATION FOR MAIN  ----------------------\n";
        for(auto it = tmain.cbegin(); it != tmain.cend(); ++it) {
            std::size_t found = it->first.find_first_of(":");
            if(found != std::string::npos) {
                logfile << "  ";
                logfile << setw(maxlen+14) << left << it->first << right << it->second << "\n";
            }
            else {
                logfile << "\n";
                logfile << setw(maxlen+16) << left << it->first << right <<  it->second << "\n";
                logfile << "-----------------------------------------------------------------------------\n";
            }
        }

        logfile << "\n\n";
        logfile << "------------------------- TIMING INFORMATION FOR THREADS  -------------------\n\n";
        logfile << "                                           Min            Max            Avg";

        auto it1 = tmin.cbegin();
        auto it2 = tmax.cbegin();
        auto it3 = tavg.cbegin();
        while(it1 != tmin.cend()) {
            std::size_t found = it1->first.find_first_of(":");
            if(found != std::string::npos) {
                logfile << "  ";
                logfile << setw(maxlen) << left << it1->first
                << setw(16) << right << it1->second
                << setw(15) << right << it2->second
                << setw(15) << right
                << it3->second/T.get_threads_per_node() << "\n";
            }
            else {
                logfile << "\n";
                logfile << setw(maxlen) << left << it1->first
                << setw(18) << right << it1->second
                << setw(15) << right << it2->second
                << setw(15) << right << it3->second/T.get_threads_per_node() 
                << "\n-----------------------------------------------------------------------------\n";
            }
            it1++;
            it2++;
            it3++;
        }

    }

}


// Temporary until C to C++ migration is completed. Use carefully!
extern "C" void *BeginRmgTimer(const char *what) {
    RmgTimer *RT = new RmgTimer(what);
    return (void *)RT;
}
extern "C" void EndRmgTimer(void *ptr) {
    RmgTimer *RT = (RmgTimer *)ptr;
    delete(RT);
}

extern "C" void CompatRmgTimerPrint(const char *outfile) {
    RmgTimer RT("Print timings");
    RT.PrintTimings(outfile);
}

