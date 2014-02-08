
#include "stdio.h"
#include "mpi.h"
#include "rmg_error.h"

void error_handler(char *message)
{

}
void rmg_error_handler(char *message)
{
    printf("%s\n", message);
    fflush (NULL);
    sleep (2);
    MPI_Abort( MPI_COMM_WORLD, 0 );
}
