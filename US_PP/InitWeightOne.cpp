/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "AtomicInterpolate.h"
#include "RmgException.h"


void InitWeightOne (SPECIES * sp, fftw_complex * rtptr, std::complex<double> *phaseptr, int ip, int l, int m, 
        fftw_plan p0_fine, fftw_plan p1_fine, fftw_plan p2_forward, fftw_plan p2_backward)
{

    int idx, idx1, size;
    double r, ax[3], bx[3];
    double t1;
    std::complex<double> *weptr, *gwptr;
    int xdim, ydim, zdim;


    double hxx = get_hxgrid() / (double) ct.nxfgrid;
    double hyy = get_hygrid() / (double) ct.nyfgrid;
    double hzz = get_hzgrid() / (double) ct.nzfgrid;
    double xside = get_xside();
    double yside = get_yside();
    double zside = get_zside();

    /* nl[xyz]fdim is the size of the non-local box in the high density grid */
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    weptr = new std::complex<double>[size];
    gwptr = new std::complex<double>[size];



    if(!weptr || !gwptr)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";



    // We assume that ion is in the center of non-local box for the localized
    // projector case. For the non-localized case it does not matter as long as
    // usage is consistent here and in GetPhase.cpp

    for(idx = 0; idx < size; idx++) weptr[idx] = 0.0;

    RmgTimer *RT0 = new RmgTimer("weight cal");
    for (int ix = 1; ix < sp->nlfdim; ix++)
    {
        double xc = (double) (ix-sp->nlfdim/2) *hxx;
        for (int iy = 1; iy < sp->nlfdim; iy++)
        {
            double yc = (double) (iy-sp->nlfdim/2) *hyy;
            for (int iz = 1; iz < sp->nlfdim; iz++)
            {
                double zc = (double) (iz-sp->nlfdim/2) *hzz;


                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);
                if(r > sp->nlradius) continue;
                to_cartesian (ax, bx);
                r += 1.0e-10;

                t1 = AtomicInterpolateInline(&sp->betalig[ip][0], r);
                idx1 = ix * sp->nlfdim * sp->nlfdim + iy * sp->nlfdim + iz;
                weptr[idx1] += Ylm(l, m, bx) * t1 ;

            }                   /* end for */

        }                       /* end for */


    }                           /* end for */

    delete RT0;
    RmgTimer *RT = new RmgTimer("weight fft_fine");
    fftw_execute_dft (p1_fine, reinterpret_cast<fftw_complex*>(weptr), reinterpret_cast<fftw_complex*>(gwptr));
    delete RT;

    PackGftoc (sp->nlfdim, sp->nlfdim, sp->nlfdim, sp->nldim, sp->nldim, sp->nldim, gwptr, weptr);

//    for(int ix = 0; ix < sp->nldim; ix++) printf("\n aaa %d %e %e", ix, weptr[ix]);

    RmgTimer *RT1 = new RmgTimer("weight fft_nldim");
    fftw_execute_dft (p2_backward, reinterpret_cast<fftw_complex*>(weptr), reinterpret_cast<fftw_complex*>(gwptr));
    delete RT1;

    RmgTimer *RT2 = new RmgTimer("weight fold");
    xdim = std::min(sp->nldim, get_NX_GRID());
    ydim = std::min(sp->nldim, get_NY_GRID());
    zdim = std::min(sp->nldim, get_NZ_GRID());
    size = xdim * ydim * zdim;

    for(idx = 0; idx < size; idx++) weptr[idx] = 0.0;
    for (int ix = 0; ix < sp->nldim; ix++)
    {
        int ixx = (ix-sp->nldim/2 + xdim/2 + 20 * xdim) % xdim;

        for (int iy = 0; iy < sp->nldim; iy++)
        {
            int iyy = (iy-sp->nldim/2 + ydim/2 + 20 * ydim) % ydim;
            for (int iz = 0; iz < sp->nldim; iz++)
            {
                int izz = (iz-sp->nldim/2 + zdim/2 + 20 * zdim) % zdim;
                idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                idx = ixx * ydim * zdim + iyy * zdim + izz;
                weptr[idx] += gwptr[idx1] * phaseptr[idx1]/double(size);
            }
        }
    }


    delete RT2;
    RmgTimer *RT3 = new RmgTimer("weight fft_forward");
    fftw_execute_dft (p2_forward, reinterpret_cast<fftw_complex*>(weptr), rtptr);
    delete RT3;

    delete []gwptr;
    delete []weptr;
}

