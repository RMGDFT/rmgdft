
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
    // Time the main thread in the zeroth slot
    static std::unordered_map<std::string, double> timings[MAX_RMG_THREADS+1];

public:
    RmgTimer(const char *what);
    ~RmgTimer(void);
    void PrintTimings(const char *outfile);

};
#else
void CompatRmgTimerPrint(const char *outfile);
void *BeginRmgTimer(const char *what);
void EndRmgTimer(void *ptr);
#endif

#endif
