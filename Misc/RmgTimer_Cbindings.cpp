#include "RmgTimer.h"
#include "transition.h"

void RmgPrintTimings(BaseGrid *G, const char *outfile, int steps);

extern "C" void CompatRmgTimerPrint(const char *outfile, int steps)
{
    RmgPrintTimings(Rmg_G, outfile, steps);
}

