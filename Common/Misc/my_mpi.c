/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"

#if MPI

void my_barrier ()
{
    MPI_Barrier (pct.img_comm);
    /*MPI_Barrier (MPI_COMM_WORLD);*/

}

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
