/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"

/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
/*This does the same as init_weight except we do it for for derivative of beta*/
/*This should be called after init_weight*/
void init_derweight (void)
{

#if !FDIFF_BETA
    int ip, prjcount, isp, size;
    SPECIES *sp;
    fftw_plan p1;
    double complex *in, *out;


    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {


        /* Get species type */
        sp = &ct.sp[isp];

        size = sp->nldim * sp->nldim * sp->nldim;

        /*This array will store forward fourier transform on the coarse grid for all betas */
        my_malloc (sp->forward_derbeta_x, 3 * sp->num_projectors * size, fftw_complex);
        if (sp->forward_derbeta_x == NULL)
            error_handler ("Could not get memory to store forward fourier transform");

        sp->forward_derbeta_y = sp->forward_derbeta_x + sp->num_projectors * size;
        sp->forward_derbeta_z = sp->forward_derbeta_y + sp->num_projectors * size;

        in = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
        out = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);

        if(!in || !out)
            error_handler ("can't allocate memory\n");

        p1 = fftw_plan_dft_3d (sp->nlfdim, sp->nlfdim, sp->nlfdim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


        prjcount = 0;
        /* Loop over radial projectors */
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            switch (sp->llbeta[ip])
            {

            case S_STATE:
                init_derweight_s (sp,
                                  &sp->forward_derbeta_x[prjcount * size],
                                  &sp->forward_derbeta_y[prjcount * size],
                                  &sp->forward_derbeta_z[prjcount * size], ip, p1);

                prjcount += 1;
                break;

            case P_STATE:
                init_derweight_p (sp,
                                  &sp->forward_derbeta_x[prjcount * size],
                                  &sp->forward_derbeta_y[prjcount * size],
                                  &sp->forward_derbeta_z[prjcount * size], ip, p1);

                prjcount += 3;
                break;

            case D_STATE:
                init_derweight_d (sp,
                                  &sp->forward_derbeta_x[prjcount * size],
                                  &sp->forward_derbeta_y[prjcount * size],
                                  &sp->forward_derbeta_z[prjcount * size], ip, p1);
                prjcount += 5;
                break;

            case F_STATE:
                error_handler ("Angular momentum state not programmed");
                break;

            default:
                error_handler ("Angular momentum state not programmed");

            }                   /* end switch */

        }                       /*end for */

        fftw_destroy_plan (p1);
        fftw_free(out);
        fftw_free(in);


    }                           /* end for */



#endif
}                               /* end get_weight */
