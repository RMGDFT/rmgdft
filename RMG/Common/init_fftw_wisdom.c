/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <fftw3.h>
#include <omp.h>
#include "main.h"
void init_fftw_wisdom (void)
{
//fftw_plan_with_nthreads(omp_get_max_threads());
#if 0
    int is;
    fftw_plan p2;
    /*fftw_plan p1; */
    for (is = 0; is < ct.num_species; is++)
    {
        /*Backward wisdom */
        p2 = fftw3d_create_plan (ct.sp[is].nldim, ct.sp[is].nldim, ct.sp[is].nldim,
                                 FFTW_BACKWARD, FFTW_MEASURE | FFTW_USE_WISDOM);

        ct.sp[is].backward_wisdom = fftw_export_wisdom_to_string ();
        if (ct.sp[is].backward_wisdom == NULL)
            error_handler ("Error when getting forward_wisdom");


        fftw_destroy_plan (p2);
        fftw_forget_wisdom ();
    }

    /*Setup this local variable to 1, to say that wisdom was succesfully allocated */
    ct.fftw_wisdom_setup = 1;
#endif
    ct.fftw_wisdom_setup = 0;
}
