/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "common_prototypes.h"
#include "main.h"

void get_weight (void)
{

    int ion, ion1, ip, coarse_size, max_size, idx, P0_BASIS;
    rmg_double_t *rtptr, *Bweight;
#if FDIFF_BETA
    rmg_double_t *rtptr_x, *rtptr_y, *rtptr_z;
    rmg_double_t *r1, *r2, *r3;
#endif
    SPECIES *sp;
    ION *iptr;
    fftw_plan p2;
    /*Pointer to the result of forward transform on the coarse grid */
    double complex *fptr;
    double complex *beptr, *gbptr;
    double complex *in, *out;

    P0_BASIS = get_P0_BASIS();

    /*maximum of nldim^3 for any species */
    max_size = ct.max_nldim * ct.max_nldim * ct.max_nldim;


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    beptr = (double complex *)fftw_malloc(sizeof(double complex) * 2 * max_size);

//    my_malloc (beptr, 2 * max_size, fftw_complex);
    if (beptr == NULL)
        error_handler ("can't allocate memory\n");

    gbptr = beptr + max_size;

#if FDIFF_BETA
    my_malloc (r1, 3 * max_size, rmg_double_t);
    r2 = r1 + max_size;
    r3 = r2 + max_size;
#endif


    for(idx = 0; idx < pct.num_tot_proj * P0_BASIS; idx++)
    {
        pct.weight[idx] = 0.0;
        pct.Bweight[idx] = 0.0;
    }
    /* Loop over ions */
    for (ion1 = 0; ion1 < pct.num_nonloc_ions; ion1++)
    {
        rtptr = &pct.weight[ion1 * ct.max_nl * P0_BASIS];
        Bweight = &pct.Bweight[ion1 * ct.max_nl * P0_BASIS];
        ion = pct.nonloc_ions_list[ion1];
        /* Generate ion pointer */
        iptr = &ct.ions[ion];


        /* Get species type */
        sp = &ct.sp[iptr->species];

        in = (double complex *)fftw_malloc(sizeof(double complex) * sp->nldim * sp->nldim * sp->nldim);
        out = (double complex *)fftw_malloc(sizeof(double complex) * sp->nldim * sp->nldim * sp->nldim);


        /*Number of grid points on which fourier transform is done (in the corse grid) */
        coarse_size = sp->nldim * sp->nldim * sp->nldim;



        //fftw_import_wisdom_from_string (sp->backward_wisdom);
        p2 = fftw_plan_dft_3d (sp->nldim, sp->nldim, sp->nldim, in, out, FFTW_BACKWARD,
                FFTW_ESTIMATE);
        //fftw_forget_wisdom ();


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
                gbptr[idx] =
                    (creal(fptr[idx]) * iptr->fftw_phase_cos[idx] + cimag(fptr[idx]) * iptr->fftw_phase_sin[idx]) +
                    (cimag(fptr[idx]) * iptr->fftw_phase_cos[idx] * 1.0I - creal(fptr[idx]) * iptr->fftw_phase_sin[idx] * 1.0I);
            }


            /*Do the backwards transform */
            fftw_execute_dft (p2, gbptr, beptr);
            /*This takes and stores the part of beta that is useful for this PE */
            assign_weight (sp, ion, beptr, rtptr, Bweight);

            /*Calculate derivative of beta */
#if FDIFF_BETA
            partial_beta_fdiff (beptr, sp->nldim, r1, r2, r3);

            assign_weight2 (sp->nldim, ion, r1, rtptr_x);
            assign_weight2 (sp->nldim, ion, r2, rtptr_y);
            assign_weight2 (sp->nldim, ion, r3, rtptr_z);
#endif



            /*Advance the temp pointers */
            fptr += coarse_size;
            rtptr += P0_BASIS;
            Bweight += P0_BASIS;
#if FDIFF_BETA
            rtptr_x += pct.idxptrlen[ion];
            rtptr_y += pct.idxptrlen[ion];
            rtptr_z += pct.idxptrlen[ion];
#endif

        }                   /*end for(ip = 0;ip < sp->num_projectors;ip++) */


        fftw_destroy_plan (p2);
        fftw_free(out);
        fftw_free(in);



    }                           /* end for */

#if FDIFF_BETA
    my_free (r1);
#endif
    fftw_free (beptr);


}                               /* end get_weight */
