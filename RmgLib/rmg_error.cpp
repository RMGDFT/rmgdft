#include "rmg_error.h"

void rmg_error_handler(const char *filename, int line, char const *message)
{
    printf("%s at LINE %d in %s.\n", message, line, filename);
    fflush (NULL);
    sleep (2);
    MPI_Abort( MPI_COMM_WORLD, 0 );
}
