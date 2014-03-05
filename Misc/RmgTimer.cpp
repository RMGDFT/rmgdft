
#include "RmgTimer.h"
#include <string>

#if 0
volatile double start_time;
volatile double end_time;
const char *sname;

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
    if(Tptr) {
       if(!Tptr->thread_timings.count(sname)) {
           std::string *tstr = new std::string(sname);
           Tptr->thread_timings[*tstr] = end_time - start_time;
       }
       else {
           std::string tstr(sname);
           Tptr->thread_timings[tstr] += end_time - start_time;
       }
    }
}
#endif

