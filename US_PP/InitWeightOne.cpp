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

void InitWeightOne (SPECIES * sp, fftw_complex * rtptr, int ip, int l, int m, fftw_plan p1)
{

    int idx, size;
    double r, ax[3], bx[3];
    double t1;
    std::complex<double> *weptr, *gwptr;

    // define functions to distiguish s, px, py, pz, ....

    double hxx = get_hxgrid() / (double) ct.nxfgrid;
    double hyy = get_hygrid() / (double) ct.nyfgrid;
    double hzz = get_hzgrid() / (double) ct.nzfgrid;
    double xside = get_xside();
    double yside = get_yside();
    double zside = get_zside();

    int xdim = sp->nlfdim;
    int ydim = sp->nlfdim;
    int zdim = sp->nlfdim;
    if(!ct.localize_projectors) {
        xdim = ct.nxfgrid * get_NX_GRID();
        ydim = ct.nxfgrid * get_NY_GRID();
        zdim = ct.nxfgrid * get_NZ_GRID();
    }

    /* nl[xyz]fdim is the size of the non-local box in the high density grid */
    size = xdim * ydim * zdim;

    weptr = new std::complex<double>[size];
    gwptr = new std::complex<double>[size];



    // We assume that ion is in the center of non-local box for the localized
    // projector case. For the non-localized case it does not matter as long as
    // usage is consistent here and in GetPhase.cpp

    for(idx = 0; idx < xdim * ydim * zdim; idx++) weptr[idx] = 0.0;

    for (int ix = 1; ix < sp->nlfdim; ix++)
    {
        int ixx = (ix-sp->nlfdim/2 + xdim/2) % xdim;
        double xc = (double) (ix-sp->nlfdim/2) *hxx;
        for (int iy = 1; iy < sp->nlfdim; iy++)
        {
            int iyy = (iy - sp->nlfdim/2 + ydim/2) % ydim;
            double yc = (double) (iy-sp->nlfdim/2) *hyy;
            for (int iz = 1; iz < sp->nlfdim; iz++)
            {
                int izz = (iz - sp->nlfdim/2 + zdim/2) % zdim;
                double zc = (double) (iz-sp->nlfdim/2) *hzz;


                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);
                if(r > sp->nlradius) continue;
                to_cartesian (ax, bx);
                r += 1.0e-10;

                t1 = AtomicInterpolateInline(&sp->betalig[ip][0], r);
                idx = ixx * ydim * zdim + iyy * zdim + izz;
                weptr[idx] += Ylm(l, m, bx) * t1;

                //if((ix*2 + nlfxdim) == 0 || (iy*2 + nlfydim) == 0 || (iz*2 + nlfzdim) == 0 ) 
                //    weptr[idx] = 0.0;

            }                   /* end for */

        }                       /* end for */


    }                           /* end for */

    fftw_execute_dft (p1, reinterpret_cast<fftw_complex*>(weptr), reinterpret_cast<fftw_complex*>(gwptr));

    pack_gftoc (sp, reinterpret_cast<fftw_complex*>(gwptr), reinterpret_cast<fftw_complex*>(rtptr));

    delete []gwptr;
    delete []weptr;
}

