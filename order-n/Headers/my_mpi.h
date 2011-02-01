/************************** SVN Revision Information **************************
 **    $Id: my_mpi.h 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#if MPI

#include <mpi.h>

/* For Linux and MPICH
 *      #include "/usr/lib/mpich/include/mpi.h"
 */

MPI_Comm COMM_KP, COMM_KPSUB1, COMM_KPSUB2;
MPI_Comm COMM_PEX, COMM_3D;

#else

int MPI_COMM_WORLD;
int COMM_KP, COMM_KPSUB1, COMM_KPSUB2;
int MPI_DOUBLE, MPI_INTEGER, MPI_SUM, MPI_MAX;

void MPI_Init (int *argc, char ***argv);

void MPI_Comm_size (int MPI_COMM_WORLD, int *mpi_nprocs);

void MPI_Comm_rank (int MPI_COMM_WORLD, int *mpi_myrank);

void MPI_Finalize ();

void MPI_Cart_create (int MPI_COMM_WORLD, int ndims, int *dims,
                      int *periods, int reorder, int *COMM_KP);

void MPI_Cart_get (int COMM_KP, int ndims, int *dims,
                   int *periods, int *coords);

void MPI_Cart_sub (int COMM_KP, int *remains, int *COMM_KPSUB1);

void MPI_Cart_rank (int COMM_KP, int *coords, int *rank);

void my_barrier ();


void MPI_Allreduce (double *in, double *out, int size, int MPI_DOUBLE,
                    int MPI_SUM, int MPI_COMM_WORLD);


#endif
