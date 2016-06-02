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

    int nlfxdim = sp->nlfdim;
    int nlfydim = sp->nlfdim;
    int nlfzdim = sp->nlfdim;
    if(!ct.localize_projectors) {
        nlfxdim = ct.nxfgrid * get_NX_GRID();
        nlfydim = ct.nxfgrid * get_NY_GRID();
        nlfzdim = ct.nxfgrid * get_NZ_GRID();
    }

    /* nl[xyz]fdim is the size of the non-local box in the high density grid */
    size = nlfxdim * nlfydim * nlfzdim;

    weptr = new std::complex<double>[size];
    gwptr = new std::complex<double>[size];


    // Next we get the radius of the projectors in terms of grid points
    int dimx = 0;
    int dimy = 0;
    int dimz = 0;
    if(!ct.localize_projectors) {
        dimx = sp->nlradius / (hxx*xside);
        dimx = 2*(dimx/2);
        dimy = sp->nlradius / (hyy*yside);
        dimy = 2*(dimy/2);
        dimz = sp->nlradius / (hzz*zside);
        dimz = 2*(dimz/2);
    }

    // We assume that ion is in the center of non-local box for the localized
    // projector case. For the non-localized case it does not matter as long as
    // usage is consistent here and in GetPhase.cpp
    int ixbegin = -nlfxdim/2;
    int ixend = ixbegin + nlfxdim;
    int iybegin = -nlfydim/2;
    int iyend = iybegin + nlfydim;
    int izbegin = -nlfzdim/2;
    int izend = izbegin + nlfzdim;
    if(!ct.localize_projectors) {
        ixbegin = -dimx/2;
        ixend = ixbegin + dimx;
        iybegin = -dimy/2;
        iyend = iybegin + dimy;
        izbegin = -dimz/2;
        izend = izbegin + dimz;
    }


    for (int ix = ixbegin; ix < ixend; ix++)
    {
        int ixx = (ix + 20 * nlfxdim) % nlfxdim;

        double xc = (double) ix *hxx;

        for (int iy = iybegin; iy < iyend; iy++)
        {
            int iyy = (iy + 20 * nlfydim) % nlfydim;
            double yc = (double) iy *hyy;

            for (int iz = izbegin; iz < izend; iz++)
            {

                int izz = (iz + 20 * nlfzdim) % nlfzdim;
                double zc = (double) iz *hzz;

                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);
                to_cartesian (ax, bx);
                r += 1.0e-10;

                t1 = AtomicInterpolateInline(&sp->betalig[ip][0], r);
                idx = ixx * nlfydim * nlfzdim + iyy * nlfzdim + izz;
                weptr[idx] = Ylm(l, m, bx) * t1;

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

