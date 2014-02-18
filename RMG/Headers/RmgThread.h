#ifndef RMG_RmgThread_H
#define RMG_RmgThread_H 1

#include "BaseThread.h"


class RmgThread : public BaseThread
{

public:
    void init_threads(void);
    void run_threads(BaseThread *s);

};

#endif


