
#ifndef RMG_Timers_H
#define RMG_Timers_H 1

#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <functional>
#include <unordered_map>
#include <string>
#include "BaseThread.h"

using namespace std;

class RmgTimer {

private:
    volatile double start_time;
    volatile double end_time;
    const char *sname;


public:
    RmgTimer(const char *fname);
    ~RmgTimer(void);

};

#endif
