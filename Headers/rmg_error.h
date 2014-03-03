#ifndef RMG_ERROR_H
#define RMG_ERROR_H 1

#if __cplusplus
#include <cstdio>
#include <unistd.h>
#include <mpi.h>

class RmgError {

public:
    void rmg_error_handler(const char *filename, int line, char const *message)
    {
        printf("%s at LINE %d in %s.\n", message, line, filename);
        fflush (NULL);
        sleep (2);
        MPI_Abort( MPI_COMM_WORLD, 0 );
    }
};

#else
void rmg_error_handler(char *message);
#endif


#endif
