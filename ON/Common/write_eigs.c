/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


  write_eigs.c


    Various functions to write out results to stdout.



*/




#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main_on.h"



/* Writes eigenvalues */
void write_eigs(STATE * states)
{
    int i;

    printf("\n  KOHN-SHAM EIGENVALUES: [eV]\n");

    for (i = 0; i < ct.num_states; i++)
        printf("  %8.3f [%4.2f]%s", states[i].eig * Ha_eV,
               states[i].occupation, ((i % 5 == 4) ? "\n" : ""));

    printf("\n\n");

}                               /* end write_eigs */
