
#include "portability.h"
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
#if (defined(_WIN32) || defined(_WIN64))
    Sleep(2000);
#else
    sleep (2);
#endif
    MPI_Abort( MPI_COMM_WORLD, 0 );
}
