/************************** SVN Revision Information **************************
 **    $Id: my_mpi.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#include "params.h"
#include "my_mpi.h"

#if MPI


void my_barrier(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
}


#endif



#if SERIAL
void MPI_Init(int *argc, char ***argv)
{
}
void MPI_Comm_size(int MPI_COMM_WORLD, int *mpi_nprocs)
{
    *mpi_nprocs = 1;
}
void MPI_Comm_rank(int MPI_COMM_WORLD, int *mpi_myrank)
{
    *mpi_myrank = 0;
}

void MPI_Finalize()
{
}

void MPI_Cart_create(int MPI_COMM_WORLD, int ndims, int *dims,
                     int *periods, int reorder, int *COMM_KP)
{
}

void MPI_Cart_get(int COMM_KP, int ndims, int *dims, int *periods, int *coords)
{
    coords[0] = 0;
    coords[1] = 0;
}

void MPI_Cart_sub(int COMM_KP, int *remains, int *COMM_KPSUB1)
{
}
void MPI_Cart_rank(int COMM_KP, int *coords, int *rank)
{
    rank = 0;
}


void my_barrier()
{
}

/*
void MPI_Allreduce(double *in, double *out, int size, int MPI_DOUBLE,
                     int MPI_SUM, int MPI_COMM_WORLD)
{
 int i;
 for(i = 0; i <size; i++) out[i] = in[i];

}
*/

#endif
