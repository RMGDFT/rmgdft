/************************** SVN Revision Information **************************
 **    $Id: my_mpi.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include "params.h"
#include "my_mpi.h"



void my_barrier (void)
{
    MPI_Barrier (MPI_COMM_WORLD);
}


