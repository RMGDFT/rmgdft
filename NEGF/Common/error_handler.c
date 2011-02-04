/************************** SVN Revision Information **************************
 **    $Id: error_handler.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#define ERROR_HANDLER_C 1

#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include "md.h"

void error_handler (int pe, char *filename, int line, char *message)
{

    printf ("\nExit from PE %d, in file %s, line %d\nError Message is: '%s'\n",
            pe, filename, line, message);
    fflush (NULL);


    char tag[10];
    sprintf( tag, "%d", ct.runflag );
    my_alloc_report( tag );

    MPI_Finalize ();
    sleep (5);
    exit (0);
}
