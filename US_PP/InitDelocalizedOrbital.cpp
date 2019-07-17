#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <complex>


#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgException.h"
#include "RmgShm.h"
#include "RmgParallelFft.h"
#include "AtomicInterpolate.h"
#include "packfuncs.h"


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void InitDelocalizedOrbital (void)
{
    RmgTimer RT0("Orbital");
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;


    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;
    double ax[3];


    // get tot number of projectors and their information
    int tot_orbitals = 0;
    
    int dimx = Rmg_G->get_PX0_GRID(1);
    int dimy = Rmg_G->get_PY0_GRID(1);
    int dimz = Rmg_G->get_PZ0_GRID(1);
    int ixstart = Rmg_G->get_PX_OFFSET(1);
    int iystart = Rmg_G->get_PY_OFFSET(1);
    int izstart = Rmg_G->get_PZ_OFFSET(1);
    int pbasis = Rmg_G->get_P0_BASIS(1);
    double vol = Rmg_L.get_omega();

    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> phase = PI * I_t;
    phase = std::exp(phase);


    std::complex<double> *weptr = new std::complex<double>[pbasis];

    RmgTimer *RT1= new RmgTimer("Orbital: phase and set");
    for (int isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        SPECIES *sp = &ct.sp[isp];

        /*Loop over all radial atomic wavefunctions to calculate num of 3D orbitals for given species */
        int orbital_count = 0;
        for (int ip = 0; ip < sp->num_atomic_waves; ip++)
        {
            if(sp->atomic_wave_oc[ip] > 0.0) {

                for(int m = 0; m < 2*sp->atomic_wave_l[ip]+1; m++)
                {
                    proj.species = isp;
                    proj.ip = ip;
                    proj.l = sp->atomic_wave_l[ip];
                    proj.m = m;
                    proj.proj_index = orbital_count;
                    proj_iter.push_back(proj);
                    orbital_count++;
                }
            }
        }

        sp->num_orbitals = orbital_count;
        tot_orbitals += orbital_count;

        /*This array will store forward fourier transform on the coarse grid for all atomic orbitals of this species */
        sp->forward_orbital = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_orbitals * pbasis * ct.num_kpts_pe);

        if (sp->forward_orbital == NULL)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    }


    double gcut = sqrt(ct.filter_factor*coarse_pwaves->gcut*tpiba2);
    RmgTimer *RT3= new RmgTimer("Orbital: proj cal");
    for(int iproj = 0; iproj < tot_orbitals; iproj++)
    {
        proj = proj_iter[iproj];
        std::complex<double> IL = std::pow(-I_t, proj.l);
        SPECIES *sp = &ct.sp[proj.species];

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            int kpt1 = kpt + pct.kstart;
            std::complex<double> *betaptr = (std::complex<double> *)&sp->forward_orbital[kpt *sp->num_orbitals * pbasis + proj.proj_index * pbasis];
            for(int idx = 0;idx < pbasis;idx++)
            {
                if(!coarse_pwaves->gmask[idx]) continue;
                weptr[idx] = std::complex<double>(0.0,0.0);
                ax[0] = coarse_pwaves->g[idx].a[0] * tpiba;
                ax[1] = coarse_pwaves->g[idx].a[1] * tpiba;
                ax[2] = coarse_pwaves->g[idx].a[2] * tpiba;

                ax[0] += ct.kp[kpt1].kvec[0];
                ax[1] += ct.kp[kpt1].kvec[1];
                ax[2] += ct.kp[kpt1].kvec[2];

                double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if(gval >= gcut) continue;
                double t1 = AtomicInterpolateInline_Ggrid(sp->atomic_wave_g[proj.ip], gval);
                weptr[idx] = IL * Ylm(proj.l, proj.m, ax) * t1;
            }

            // Shift atom to the center instead of corner.
            for(int ix = 0; ix < dimx; ix++)
            {
                for(int iy = 0; iy < dimy; iy++)
                {
                    for(int iz = 0; iz < dimz; iz++)
                    {
                        int idx = ix * dimy * dimz + iy * dimz + iz;
                        weptr[idx] *= std::pow(phase, ix + ixstart + iy + iystart + iz + izstart);
                    }
                }
            }
            for(int idx=0;idx < pbasis;idx++) betaptr[idx] =  weptr[idx]/vol;

        }

    }  // end for
    delete RT3;

    delete RT1;
    delete [] weptr;

} /* end InitDelocalizedOrbital */


