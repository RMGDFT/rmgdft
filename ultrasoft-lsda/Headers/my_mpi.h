/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#if MPI

#include <mpi.h>

#else
void MPI_Abort ();
void MPI_Finalize ();
void my_barrier ();
#endif
