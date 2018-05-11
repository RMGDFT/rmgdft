
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include "rmg_error.h"
#include "rmg_mangling.h"
#include "portability.h"
#include "params.h"
#include "Scalapack.h"
#include "pe_control.h"

#define         errore          RMG_FC_GLOBAL(errore, ERRORE)


void error_handler(char *message)
{

}
void rmg_error_handler(char *message)
{
    if(pct.gridpe == 0) printf("%s\n", message);
    fflush (NULL);
#if (defined(_WIN32) || defined(_WIN64))
    Sleep(2000);
#else
    sleep (2);
#endif

}

// Here for fortran routines.
void errore(char *where, char *message, int ierr, int where_len, int message_len)
{
  char tbuf[1000];
  memset(tbuf, 0, sizeof(tbuf));

  if(((where_len + message_len) > sizeof(tbuf)) || (where_len < 0) || (message_len < 0)) {
     printf("Unknown issue printing error message from fortran routines\n"); 
     raise(SIGTERM);
  }

  strncpy(tbuf, message, message_len);
  strncpy(&tbuf[message_len], " in ", 4);
  strncpy(&tbuf[message_len + 4], where, where_len);

  printf("%s\n", tbuf);
  fflush (NULL);
#if (defined(_WIN32) || defined(_WIN64))
    Sleep(2000);
#else
    sleep (2);
#endif
    MPI_Abort( MPI_COMM_WORLD, 0 );

}
