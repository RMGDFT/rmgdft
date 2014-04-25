/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"

void get_derweight (int ion, rmg_double_t * beta_x, rmg_double_t * beta_y, rmg_double_t * beta_z, ION * iptr,
                    fftw_plan p2)
{

#if !FDIFF_BETA
    int ip, coarse_size, idx, nh;
    SPECIES *sp;
    rmg_double_t *rtptr_x, *rtptr_y, *rtptr_z;
    /*Pointer to the result of forward transform on the coarse grid */
    complex *fptr_x, *fptr_y, *fptr_z;
    complex *beptr, *gbptr;



    /* Get species type */
    sp = &ct.sp[iptr->species];

    /*Number of grid points on which fourier transform is done (in the corse grid) */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;

    nh = pct.idxptrlen[ion];


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    my_malloc (beptr, 2 * coarse_size, fftw_complex);
    if (beptr == NULL)
        error_handler ("can't allocate memory\n");

    gbptr = beptr + coarse_size;



    /*Temporary pointer to the already calculated forward transform */
    fptr_x = sp->forward_derbeta_x;
    fptr_y = sp->forward_derbeta_y;
    fptr_z = sp->forward_derbeta_z;

    rtptr_x = beta_x;
    rtptr_y = beta_y;
    rtptr_z = beta_z;


    /* Loop over radial projectors */
    for (ip = 0; ip < sp->num_projectors; ip++)
    {


    /*************** X ********************/
        /*Apply the phase factor */
        for (idx = 0; idx < coarse_size; idx++)
        {
            gbptr[idx] =
                (creal(fptr_x[idx]) * iptr->fftw_phase_cos[idx] +
                cimag(fptr_x[idx]) * iptr->fftw_phase_sin[idx]) +
                (cimag(fptr_x[idx]) * iptr->fftw_phase_cos[idx] -
                creal(fptr_x[idx]) * iptr->fftw_phase_sin[idx]) * 1.0I;
        }

        /*Do the backwards transform */
        fftw_execute_dft (p2, gbptr, beptr);
        /*This takes and stores the part of beta that is useful for this PE */
        assign_derweight (sp, ion, beptr, rtptr_x);


    /*************** Y ********************/
        /*Apply the phase factor */
        for (idx = 0; idx < coarse_size; idx++)
        {
            gbptr[idx] =
                (creal(fptr_y[idx]) * iptr->fftw_phase_cos[idx] +
                cimag(fptr_y[idx]) * iptr->fftw_phase_sin[idx]) +
                (cimag(fptr_y[idx]) * iptr->fftw_phase_cos[idx] -
                creal(fptr_y[idx]) * iptr->fftw_phase_sin[idx]) * 1.0I;
        }

        /*Do the backwards transform */
        fftw_execute_dft (p2, gbptr, beptr);
        /*This takes and stores the part of beta that is useful for this PE */
        assign_derweight (sp, ion, beptr, rtptr_y);


    /*************** Z ********************/
        /*Apply the phase factor */
        for (idx = 0; idx < coarse_size; idx++)
        {
            gbptr[idx] =
                (creal(fptr_z[idx]) * iptr->fftw_phase_cos[idx] +
                cimag(fptr_z[idx]) * iptr->fftw_phase_sin[idx]) +
                (cimag(fptr_z[idx]) * iptr->fftw_phase_cos[idx] -
                creal(fptr_z[idx]) * iptr->fftw_phase_sin[idx]) * 1.0I;
        }

        /*Do the backwards transform */
        fftw_execute_dft (p2, gbptr, beptr);
        /*This takes and stores the part of beta that is useful for this PE */
        assign_derweight (sp, ion, beptr, rtptr_z);



        /*Advance the temp pointers */
        fptr_x += coarse_size;
        fptr_y += coarse_size;
        fptr_z += coarse_size;
        rtptr_x += nh;
        rtptr_y += nh;
        rtptr_z += nh;

    }                           /*end for(ip = 0;ip < sp->num_projectors;ip++) */





    my_free (beptr);


#endif
}                               /* end get_weight */
