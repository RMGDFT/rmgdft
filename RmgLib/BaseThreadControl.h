#ifndef RMG_BaseThreadControl_H
#define RMG_BaseThreadControl_H 1

class BaseThreadControl {

public:

    // Base tag
    int basetag;

    // Thread ID number assigned by us
    int tid;

    // Pointer to project specific data structure
    void *pptr;

};


#endif
