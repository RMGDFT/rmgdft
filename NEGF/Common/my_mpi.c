/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include "params.h"
#include "my_mpi.h"



void my_barrier (void)
{
    MPI_Barrier (MPI_COMM_WORLD);
}


