#ifndef PORTABILITY_H_INCLUDED
#define PORTABILITY_H

#if (defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64))
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #include <windows.h>
    #define snprintf _snprintf
    typedef int mode_t;
    #include <io.h>
#else
    #include <unistd.h>
#endif
#endif
