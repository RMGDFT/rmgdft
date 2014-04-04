#include "rmg_error.h"

static void (*rmgerrfunc)(const char *filename, int line, char const *message) = NULL;

void RmgRegisterErrorHandler(void (*func)(const char *filename, int line, char const *message))
{
    rmgerrfunc = func;
}

void rmg_error_handler(const char *filename, int line, char const *message)
{
    if(!rmgerrfunc) {
        printf("%s at LINE %d in %s.\n", message, line, filename);
        fflush (NULL);
        sleep (2);
        MPI_Abort( MPI_COMM_WORLD, 0 );
    }
    rmgerrfunc(filename, line, message);
}
