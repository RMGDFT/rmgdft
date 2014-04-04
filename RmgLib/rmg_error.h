#ifndef RMG_ERROR_H
#define RMG_ERROR_H 1

#if __cplusplus
#include <cstdio>
#include <unistd.h>
#include <mpi.h>

void rmg_error_handler(const char *filename, int line, char const *message);
void RmgRegisterErrorHandler(void (*func)(const char *filename, int line, char const *message));

#else
void rmg_error_handler(char *message);
#endif


#endif
