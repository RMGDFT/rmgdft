/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include "Atomic.h"

#include "AtomicInterpolate.h"

void InitLocalForward ()
{
    int isp, size;
    SPECIES *sp;
    fftw_plan p1;
    fftw_complex *in, *out;
    std::complex<double> *vnuc_ptr, *rhoc_ptr, *rhocore_ptr, *gwptr;
    int idx, ix, iy, iz, ibegin, iend;
    double r, ax[3], xc, yc, zc;
    double t1, hxx, hyy, hzz;
    double Zv, rc, rc2, rcnorm;
    Atomic *A = new Atomic();

    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {


        /* Get species type */
        sp = &ct.sp[isp];
        Zv = sp->zvalence;
        rc = sp->rc;
        rc2 = rc * rc;
        rcnorm = rc * rc * rc * pow (PI, 1.5);
        rcnorm = 1.0 / rcnorm;

        size = sp->ldim * sp->ldim * sp->ldim;

        /*This array will store forward fourier transform on the coarse grid for all betas */
        sp->forward_vnuc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);
        sp->forward_rhoc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);
        sp->forward_rhocore = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);

        if (sp->forward_beta == NULL)
        {
            printf("\n Could not get memory to store forward fourier transform, vnuc");
            fflush(NULL);
            exit(0);
        }


        //  creat plan on fine grid 

        size = sp->ldim_fine * sp->ldim_fine * sp->ldim_fine;

        vnuc_ptr = new std::complex<double>[size];
        rhoc_ptr = new std::complex<double>[size];
        rhocore_ptr = new std::complex<double>[size];
        gwptr = new std::complex<double>[size];


        in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);
        out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);

        if(!in || !out)
        {
            printf("\n Could not get memory for in and put in InitLocalForward");
            fflush(NULL);
            exit(0);
        }

        p1 = fftw_plan_dft_3d (sp->ldim_fine, sp->ldim_fine, sp->ldim_fine, in, out, FFTW_FORWARD, FFTW_MEASURE);



        hxx = get_hxgrid() / (double) ct.nxfgrid;
        hyy = get_hygrid() / (double) ct.nyfgrid;
        hzz = get_hzgrid() / (double) ct.nzfgrid;

        /*We assume that ion is in the center of non-local box */
        //    ibegin = -(sp->nldim / 2) * ct.nxfgrid;

        ibegin = -sp->ldim_fine/2;
        iend = ibegin + sp->ldim_fine;
        int ixx, iyy, izz;
        for (ix = ibegin; ix < iend; ix++)
        {
            ixx = ix;
            if (ixx < 0) ixx = ix + sp->ldim_fine;
            xc = (double) ix *hxx;

            for (iy = ibegin; iy < iend; iy++)
            {
                iyy = iy;
                if (iyy < 0) iyy = iy + sp->ldim_fine;
                yc = (double) iy *hyy;

                for (iz = ibegin; iz < iend; iz++)
                {

                    izz = iz;
                    if (izz < 0) izz = iz + sp->ldim_fine;
                    zc = (double) iz *hzz;

                    ax[0] = xc;
                    ax[1] = yc;
                    ax[2] = zc;

                    r = metric (ax);

                    idx = ixx *sp->ldim_fine * sp->ldim_fine + iyy * sp->ldim_fine + izz;


                    t1= AtomicInterpolateInline (&sp->localig[0], r);
                    vnuc_ptr[idx] = t1;

                    t1 = Zv * exp (-r * r / rc2) * rcnorm;
                    rhoc_ptr[idx] = t1;


                    if (sp->nlccflag)
                    {

                        t1 = AtomicInterpolateInline (&sp->rhocorelig[0], r);
                        rhocore_ptr[idx] = t1;

                    }

                }                           /* end for */
            }
        }



        fftw_execute_dft (p1, reinterpret_cast<fftw_complex*>(vnuc_ptr), reinterpret_cast<fftw_complex*>(gwptr));
        A->PackFine2Rhogrid(gwptr, sp->ldim_fine, (std::complex<double>*)sp->forward_vnuc, sp->ldim); 

        fftw_execute_dft (p1, reinterpret_cast<fftw_complex*>(rhoc_ptr), reinterpret_cast<fftw_complex*>(gwptr));
        A->PackFine2Rhogrid(gwptr, sp->ldim_fine, (std::complex<double> *)sp->forward_rhoc, sp->ldim); 

        if (sp->nlccflag)
        {
            fftw_execute_dft (p1, reinterpret_cast<fftw_complex*>(rhocore_ptr), reinterpret_cast<fftw_complex*>(gwptr));
            A->PackFine2Rhogrid(gwptr, sp->ldim_fine, (std::complex<double> *)sp->forward_rhocore, sp->ldim); 

        }

        fftw_destroy_plan (p1);
        fftw_free(out);
        fftw_free(in);

        delete [] vnuc_ptr;
        delete [] rhoc_ptr;
        delete [] rhocore_ptr;
        delete [] gwptr;



    }                           /* end for */
}

