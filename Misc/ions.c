/*

  Data and functions related to ions.

*/

#include "common_prototypes.h"
#include "typedefs.h"


ION *get_ion(int ion)
{
    return ct.ions[ion];
}


