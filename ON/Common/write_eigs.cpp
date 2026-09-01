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
#include "main.h"
#include "transition.h"
#include "prototypes_on.h"



/* Writes eigenvalues */
void write_eigs(STATE * states, double *kpt)
{
    int i;

    rmg_printf("\n  KOHN-SHAM EIGENVALUES: [eV] at K =[%f %f %f]\n", kpt[0]/(2.0*PI), kpt[1]/(2.0*PI),kpt[2]/(2.0*PI));

    for (i = 0; i < ct.num_states; i++)
        rmg_printf("  %8.3f [%4.2f]%s", states[i].eig[0] * Ha_eV,
                states[i].occupation[0], ((i % 5 == 4) ? "\n" : ""));
    if(ct.spin_flag)
    {
        rmg_printf("\n Eigs for oppsite spin\n");
        for (i = 0; i < ct.num_states; i++)
            rmg_printf("  %8.3f [%4.2f]%s", states[i].eig[1] * Ha_eV,
                    states[i].occupation[1], ((i % 5 == 4) ? "\n" : ""));
    }

    rmg_printf("\n\n");

}                               /* end write_eigs */
