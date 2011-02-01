/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"

#if MPI

#  if (AIX || LINUX || IRIX || XT3)
void my_barrier ()
{
    MPI_Barrier (pct.grid_comm);
    /*MPI_Barrier (MPI_COMM_WORLD);*/

}
#  endif

#else

void MPI_Finalize ()
{
}

void my_barrier ()
{
}

/*
void MPI_Allreduce(double *in, double *out, int size, int MPI_DOUBLE,
                     int MPI_SUM, int pct.grid_comm)
{
 int i;
 for(i = 0; i <size; i++) out[i] = in[i];

}
*/

#endif
