
#ifndef RMG_Timers_H
#define RMG_Timers_H 1



#if __cplusplus

#include <unordered_map>
#include <fstream>
#include "BaseThread.h"
#include "BaseGrid.h"


class RmgTimer {

private:
    volatile double start_time;
    volatile double end_time;
    const char *sname;

public:
    RmgTimer(const char *what);
    ~RmgTimer(void);
    std::unordered_map<std::string, double> *get_map(void);

    // The main thread of execution is timed in the zeroth slot. Data structure public
    // is public so that the main program can implement it's own output routines.
    static std::unordered_map<std::string, double> timings[MAX_RMG_THREADS+1];

};
#else
void *BeginRmgTimer(const char *what);
void EndRmgTimer(void *ptr);
#endif

#endif
