/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void get_weight (void)
{

    int ion, ion1, ip, coarse_size, max_size, idx;
    REAL *rtptr;
#if FDIFF_BETA
    REAL *rtptr_x, *rtptr_y, *rtptr_z;
    REAL *r1, *r2, *r3;
#endif
    SPECIES *sp;
    ION *iptr;
    fftwnd_plan p2;
    /*Pointer to the result of forward transform on the coarse grid */
    fftw_complex *fptr;
    fftw_complex *beptr, *gbptr;


    /*maximum of nldim^3 for any species */
    max_size = ct.max_nldim * ct.max_nldim * ct.max_nldim;


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    my_malloc (beptr, 2 * max_size, fftw_complex);
    if (beptr == NULL)
        error_handler ("can't allocate memory\n");

    gbptr = beptr + max_size;

#if FDIFF_BETA
    my_malloc (r1, 3 * max_size, REAL);
    r2 = r1 + max_size;
    r3 = r2 + max_size;
#endif


     rtptr = pct.weight;
    /* Loop over ions */
    for (ion1 = 0; ion1 < pct.num_nonloc_ions; ion1++)
    {
        ion = pct.nonloc_ions_list[ion1];
        /* Generate ion pointer */
        iptr = &ct.ions[ion];


        /* Get species type */
        sp = &ct.sp[iptr->species];

        /*Number of grid points on which fourier transform is done (in the corse grid) */
        coarse_size = sp->nldim * sp->nldim * sp->nldim;



        fftw_import_wisdom_from_string (sp->backward_wisdom);
        p2 = fftw3d_create_plan (sp->nldim, sp->nldim, sp->nldim, FFTW_BACKWARD,
                FFTW_USE_WISDOM);
        fftw_forget_wisdom ();


        /*Calculate the phase factor */
        find_phase (sp->nldim, iptr->nlcrds, iptr->fftw_phase_sin, iptr->fftw_phase_cos);


        /*Temporary pointer to the already calculated forward transform */
        fptr = sp->forward_beta;

        /*Pointer to where calculated beta will be stored */
#if FDIFF_BETA
        rtptr_x = pct.weight_derx[ion];
        rtptr_y = pct.weight_dery[ion];
        rtptr_z = pct.weight_derz[ion];
#endif


        /* Loop over radial projectors */
        for (ip = 0; ip < sp->num_projectors; ip++)
        {


            /*Apply the phase factor */
            for (idx = 0; idx < coarse_size; idx++)
            {
                gbptr[idx].re =
                    fptr[idx].re * iptr->fftw_phase_cos[idx] +
                    fptr[idx].im * iptr->fftw_phase_sin[idx];
                gbptr[idx].im =
                    fptr[idx].im * iptr->fftw_phase_cos[idx] -
                    fptr[idx].re * iptr->fftw_phase_sin[idx];
            }


            /*Do the backwards transform */
            fftwnd_one (p2, gbptr, beptr);
            /*This takes and stores the part of beta that is useful for this PE */
            assign_weight (sp, ion, beptr, rtptr);

            /*Calculate derivative of beta */
#if FDIFF_BETA
            partial_beta_fdiff (beptr, sp->nldim, r1, r2, r3);

            assign_weight2 (sp->nldim, ion, r1, rtptr_x);
            assign_weight2 (sp->nldim, ion, r2, rtptr_y);
            assign_weight2 (sp->nldim, ion, r3, rtptr_z);
#endif



            /*Advance the temp pointers */
            fptr += coarse_size;
            rtptr += pct.P0_BASIS;
#if FDIFF_BETA
            rtptr_x += pct.idxptrlen[ion];
            rtptr_y += pct.idxptrlen[ion];
            rtptr_z += pct.idxptrlen[ion];
#endif

        }                   /*end for(ip = 0;ip < sp->num_projectors;ip++) */


        fftwnd_destroy_plan (p2);


    }                           /* end for */

#if FDIFF_BETA
    my_free (r1);
#endif
    my_free (beptr);


}                               /* end get_weight */
