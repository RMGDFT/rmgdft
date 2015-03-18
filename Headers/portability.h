#ifndef PORTABILITY_H_INCLUDED
#define PORTABILITY_H

#if (defined(_WIN32) || defined(_WIN64))
    #define snprintf _snprintf
#endif

#endif
