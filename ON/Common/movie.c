/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "prototypes_on.h"

void movie(FILE * mfp)
{

    int ion;
    ION *iptr;

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];
        fprintf(mfp, " %d %f %f %f %14.10f %14.10f %14.10f\n", iptr->species,
                iptr->crds[0],
                iptr->crds[1],
                iptr->crds[2], iptr->velocity[0], iptr->velocity[1], iptr->velocity[2]);

    }
    fprintf(mfp, "/end\n");

}                               /* end of movie */
