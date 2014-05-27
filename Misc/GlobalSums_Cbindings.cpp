
#include "GlobalSums.h"
#include "common_prototypes.h"


extern "C" void init_global_sums(void) 
{
   GlobalSumsInit();
}

extern "C" void global_sums (rmg_double_t * vect, int *length, MPI_Comm comm)
{
    GlobalSums<double>(vect, *length, comm);
}
