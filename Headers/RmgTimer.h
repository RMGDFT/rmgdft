
#ifndef RMG_Timers_H
#define RMG_Timers_H 1



#if __cplusplus

#include <sys/time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <functional>
#include <unordered_map>
#include "BaseThread.h"
using namespace std;

class RmgTimer {

private:
    volatile double start_time;
    volatile double end_time;
    const char *sname;
    static std::unordered_map<std::string, double> timings[MAX_RMG_THREADS+1];  // Time the main thread in the last slot


public:
    RmgTimer(const char *fname);
    ~RmgTimer(void);
    void PrintTimings(void);

};
#else
void CompatRmgTimerPrint(void);
#endif

#endif
