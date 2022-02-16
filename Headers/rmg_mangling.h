#ifndef RMG_mangling_in_H
#define RMG_mangling_in_H 1

#include "rmg_mangling_in.h"

// A bit of a hack but cmake can't seem to get this right for cray fortran
#if RMG_CRAY_FORTRAN
    /* Mangling for Fortran module symbols without underscores. */
    #define RMG_FC_MODULE(mod_name,name, mod_NAME,NAME) name##$##mod_name##_

    /* Mangling for Fortran module symbols with underscores. */
    #define RMG_FC_MODULE_(mod_name,name, mod_NAME,NAME) name##$##mod_name##_
#endif

#endif

