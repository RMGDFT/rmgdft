/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"
#include "RmgShm.h"

/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void init_weight (void)
{

    int ip, prjcount, isp, size;
    SPECIES *sp;
    fftw_plan p1;
    double complex *in, *out;
    int maxl = 0;

    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {


        /* Get species type */
        sp = &ct.sp[isp];


        /*Loop over all betas to calculate num of projectors for given species */
        prjcount = 0;
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            switch (sp->llbeta[ip])
            {

            case S_STATE:
                prjcount += 1;
                break;

            case P_STATE:
                prjcount += 3;
                if(maxl < 1) maxl = 1;
                break;

            case D_STATE:
                prjcount += 5;
                if(maxl < 2) maxl = 2;
                break;

            case F_STATE:
                error_handler ("Angular momentum state not programmed");
                break;

            default:
                error_handler ("Angular momentum state not programmed");

            }                   /* end switch */

        }                       /*end for(ip = 0;ip < sp->nbeta;ip++) */

        /*Store number of projectors for given species */
        sp->num_projectors = prjcount;



        size = sp->nldim * sp->nldim * sp->nldim;

        /*This array will store forward fourier transform on the coarse grid for all betas */
        bool use_shared = false;
        if(pct.procs_per_host > (2*maxl+1)) {
            char sname[256];
            snprintf(sname, sizeof(sname), "RMG_ForwardBeta_%s", sp->atomic_symbol);
            sp->forward_beta = (double complex *)AllocSharedMemorySegment(sname, sizeof(double complex) * sp->num_projectors * size);
            if(sp->forward_beta) use_shared = true;
        }
        if(!sp->forward_beta) {
            sp->forward_beta = (double complex *)fftw_malloc(sizeof(double complex) * sp->num_projectors * size);
        }

        if (sp->forward_beta == NULL)
            error_handler ("Could not get memory to store forward fourier transform");



        /*This is something we need to do only once per species, so do not use wisdom */
        in = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
        out = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);

        if(!in || !out)
            error_handler ("can't allocate memory\n");
        p1 = fftw_plan_dft_3d (sp->nlfdim, sp->nlfdim, sp->nlfdim, in, out, FFTW_FORWARD, FFTW_MEASURE);
        
        prjcount = 0;
        /* Loop over radial projectors */
        for (ip = 0; ip < sp->nbeta; ip++)
        {
            switch (sp->llbeta[ip])
            {

            case S_STATE:
                init_weight_s (sp, &sp->forward_beta[prjcount * size], ip, p1, use_shared);

                prjcount += 1;
                break;

            case P_STATE:
                init_weight_p (sp, &sp->forward_beta[prjcount * size], ip, p1, use_shared);

                prjcount += 3;
                break;

            case D_STATE:
                init_weight_d (sp, &sp->forward_beta[prjcount * size], ip, p1, use_shared);
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

}                               /* end get_weight */
