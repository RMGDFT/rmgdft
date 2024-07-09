
// Function prototype is included in Gpufuncs.h which needs to be minimal because
// hip and cuda compilers have trouble with some of the general header defs.
#include "Gpufuncs.h"
#include "BaseThread.h"
#include "main.h"

int getThreadId(void)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    return tid;
}

int getThreadNum(void)
{
    return std::max(ct.OMP_THREADS_PER_NODE, ct.MG_THREADS_PER_NODE);
}

