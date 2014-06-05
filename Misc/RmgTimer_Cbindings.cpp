#include <sys/time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <functional>
#include <string>
#include <iomanip>
#include <algorithm>
#include "RmgTimer.h"

#include "main.h"
#include "transition.h"

void RmgPrintTimings(BaseGrid *G, const char *outfile, int steps);

extern "C" void CompatRmgTimerPrint(const char *outfile, int steps)
{
    RmgPrintTimings(Rmg_G, outfile, steps);
}

