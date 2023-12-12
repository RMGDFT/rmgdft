#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "AtomicInterpolate.h"
#include "RmgException.h"


void InitWeightOne (SPECIES * sp, fftw_complex * rtptr, std::complex<double> *phaseptr, int ip, int l, int m)
{

    double ax[3];
    double t1;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double scale = 2.0 * PI / sp->prj_pwave->L->celldm[0];

    /* nl[xyz]fdim is the size of the non-local box in the high density grid */
    int size = sp->nldim * sp->nldim * sp->nldim;
    size = std::max(size, get_NX_GRID() * get_NY_GRID() * get_NZ_GRID());

    std::complex<double> *weptr = new std::complex<double>[size]();
    std::complex<double> *gwptr = new std::complex<double>[size];


    if(!weptr || !gwptr)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    int ixx, iyy, izz;
    double gval, gcut;

    int xdim = sp->nldim;
    int ydim = sp->nldim;
    int zdim = sp->nldim;

    gcut = sqrt(coarse_pwaves->gcut*tpiba2); // pwave structures store the squared magnitude

    double vol = Rmg_L.get_omega() * xdim * ydim * zdim / (get_NX_GRID() * get_NY_GRID() * get_NZ_GRID()); 

    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> IL = std::pow(-I_t, l);

//  shift the atom to center with a phase 2pi *(xdim/2)/xdim

    std::complex<double> phase;
    
    phase = 2.0 * PI * ((xdim+1)/2)/xdim * I_t;
    phase = std::exp(phase);

    for (int ix = -xdim/2; ix < xdim/2+1; ix++)
    {
        ixx = (ix + xdim)%xdim;
        for (int iy = -ydim/2; iy < ydim/2+1; iy++)
        {
            iyy = (iy + ydim)%ydim;
            for (int iz = -zdim/2; iz < zdim/2+1; iz++)
            {
                izz = (iz + zdim)%zdim;
                int idx1 = ixx * ydim * zdim + iyy * zdim + izz;

                ax[0] = scale * sp->prj_pwave->g[idx1].a[0];
                ax[1] = scale * sp->prj_pwave->g[idx1].a[1];
                ax[2] = scale * sp->prj_pwave->g[idx1].a[2];
                gval = scale*sqrt(sp->prj_pwave->gmags[idx1]);

                if(gval >= gcut) continue;
                t1 = AtomicInterpolateInline_Ggrid(&sp->beta_g[ip][0], gval);

                weptr[idx1] += IL * Ylm(l, m, ax) * t1/vol;
            }
        }
    }

// shift atom to the center instead of corner
    for(int ix = 0; ix < xdim; ix++)
    for(int iy = 0; iy < ydim; iy++)
    for(int iz = 0; iz < zdim; iz++)
    {
        int idx = ix * ydim * zdim + iy * zdim + iz;
        weptr[idx] *= std::pow(phase, ix+iy+iz);
//        if( (ix + iy + iz) %2 ) weptr[idx] *=-1.0;
    }



    //    for(int ix = 0; ix < sp->nldim; ix++) printf("\n aaa %d %e %e", ix, weptr[ix]);

    RmgTimer *RT1 = new RmgTimer("weight fft_nldim");
    sp->prj_pwave->FftInverse(weptr, gwptr);
    delete RT1;

    RmgTimer *RT2 = new RmgTimer("weight fold");
    if(!ct.localize_projectors)
    {
        xdim = get_NX_GRID();
        ydim = get_NY_GRID();
        zdim = get_NZ_GRID();
    }

    size = xdim * ydim * zdim;

    for(int idx = 0; idx < size; idx++) weptr[idx] = 0.0;
    for (int ix = 0; ix < sp->nldim; ix++)
    {
        int ixx = (ix-sp->nldim/2 + xdim/2 + 20 * xdim) % xdim;

        for (int iy = 0; iy < sp->nldim; iy++)
        {
            int iyy = (iy-sp->nldim/2 + ydim/2 + 20 * ydim) % ydim;
            for (int iz = 0; iz < sp->nldim; iz++)
            {
                int izz = (iz-sp->nldim/2 + zdim/2 + 20 * zdim) % zdim;
                int idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                int idx = ixx * ydim * zdim + iyy * zdim + izz;
                weptr[idx] += gwptr[idx1] * phaseptr[idx1]/double(size);
            }
        }
    }


    delete RT2;
    RmgTimer *RT3 = new RmgTimer("weight fft_forward");
    sp->prj_pwave->FftForward(weptr, (std::complex<double> *)rtptr);
    delete RT3;

    // shift atom to the center instead of corner

    delete []gwptr;
    delete []weptr;
}

