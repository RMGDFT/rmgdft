/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

/*Releases memory and forgets fftw wisdom*/
void destroy_fftw_wisdom (void)
{

    int is;

    /*Release fftw memory, but only if memory was previously allocated
     * ct.fftw_widom_setup measures whether wisdom was set or not, 
     * Iinit_fftw_wisdom sets this variable to 1, otherwise it is 0*/
    if (ct.fftw_wisdom_setup)
    {

        for (is = 0; is < ct.num_species; is++)
        {
            if (ct.sp[is].backward_wisdom != NULL)
                fftw_free (ct.sp[is].backward_wisdom);
            else
                error_handler ("ct.[is].backward_wisdom should be set to something");
        }
    }

    else
        printf ("\n\n PE:%d: NOT releasing fftw wisdom, it was NOT setup before", pct.gridpe);

}

 /*EOF*/
