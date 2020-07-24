#include <float.h>
#include <sys/stat.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <complex>

#include <fcntl.h>


#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgException.h"
#include "RmgParallelFft.h"
#include "AtomicInterpolate.h"
#include "packfuncs.h"


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void InitDelocalizedWeight_onek (int kindex, double kvec[3], Pw &pwave)
{
    RmgTimer RT0("Weight");


    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> *IL =  new std::complex<double>[ct.max_l+2];
    for(int L = 0; L <  ct.max_l+2; L++) IL[L] = std::pow(-I_t, L);
    std::complex<double> phase = PI * I_t;
    phase = std::exp(phase);

    int lmax = ct.max_l + 1;
    int num_lm = (lmax + 1) * (lmax + 1);
    int num_LM2 = (2*lmax + 1) * (2*lmax + 1);

    int *lpx = new int[num_lm * num_lm];
    int *lpl = new int[num_lm * num_lm  * num_LM2];
    double *ap = new double[num_LM2 * num_lm * num_lm];
    //ylm = new double[25];
    mkdir ("PROJECTORS", S_IRWXU);

    InitClebschGordan(lmax, ap, lpx, lpl);

    int pbasis = pwave.Grid->get_P0_BASIS(1);
    std::complex<double> *forward_beta;
    std::complex<double> *forward_beta_r[3];
    forward_beta = new std::complex<double>[ct.max_nl * pbasis];
    forward_beta_r[0] = new std::complex<double>[ct.max_nl * pbasis];
    forward_beta_r[1] = new std::complex<double>[ct.max_nl * pbasis];
    forward_beta_r[2] = new std::complex<double>[ct.max_nl * pbasis];

    for (int isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        SPECIES *sp = &Species[isp];


        for(int iproj = 0; iproj < sp->nh; iproj++)
        {
            int dimx = pwave.Grid->get_PX0_GRID(1);
            int dimy = pwave.Grid->get_PY0_GRID(1);
            int dimz = pwave.Grid->get_PZ0_GRID(1);
            int ixstart = pwave.Grid->get_PX_OFFSET(1);
            int iystart = pwave.Grid->get_PY_OFFSET(1);
            int izstart = pwave.Grid->get_PZ_OFFSET(1);
            double vol = pwave.L->get_omega();
            double tpiba = 2.0 * PI / Rmg_L.celldm[0];
            double tpiba2 = tpiba * tpiba;
            double gcut = sqrt(ct.filter_factor*pwave.gcut*tpiba2);
            double ax[3];

            int proj_l = sp->nhtol[iproj];
            int proj_m = sp->nhtom[iproj];
            int proj_ip = sp->indv[iproj];

            size_t index_ptr = iproj * pbasis;
            std::complex<double> *betaptr = (std::complex<double> *)&forward_beta[index_ptr];
            std::complex<double> *betaptr_r[3];

            betaptr_r[0] = (std::complex<double> *)&forward_beta_r[0][index_ptr];
            betaptr_r[1] = (std::complex<double> *)&forward_beta_r[1][index_ptr];
            betaptr_r[2] = (std::complex<double> *)&forward_beta_r[2][index_ptr];

            for(int idx = 0;idx < pbasis;idx++)
            {
                betaptr[idx] = 0.0;
                betaptr_r[0][idx] = 0.0;
                betaptr_r[1][idx] = 0.0;
                betaptr_r[2][idx] = 0.0;
            }

            for(int idx = 0;idx < pbasis;idx++)
            {
                if(!pwave.gmask[idx]) continue;
                ax[0] = pwave.g[idx].a[0] * tpiba;
                ax[1] = pwave.g[idx].a[1] * tpiba;
                ax[2] = pwave.g[idx].a[2] * tpiba;

                ax[0] += kvec[0];
                ax[1] += kvec[1];
                ax[2] += kvec[2];

                double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if(gval >= gcut) continue;
                double t1 = AtomicInterpolateInline_Ggrid(sp->beta_g[proj_ip], gval);
                betaptr[idx] = IL[proj_l] * Ylm(proj_l, proj_m, ax) * t1;

                // l2m_i: l*l + m for the first angular momentum
                int l2mi = proj_l * proj_l + proj_m;
                for(int l2mj = 1; l2mj < 4; l2mj++)     // index for cubic harmonics x, y, z
                {
                    for (int LM = 0; LM < lpx[l2mi * num_lm + l2mj]; LM++)
                    {
                        int L2M = lpl[(l2mi * num_lm + l2mj) * num_LM2 + LM];   // L*L + M for one LM harmonic function 

                        int L, M;
                        if(L2M == 0)
                            L = 0;
                        else if (L2M < 4)
                            L = 1;
                        else if (L2M < 9)
                            L = 2;
                        else
                            L = (int)sqrt(L2M + 0.1);

                        M = L2M - L * L;

                        double t2 = AtomicInterpolateInline_Ggrid(sp->rbeta_g[proj_ip][L], gval);
                        betaptr_r[l2mj-1][idx] += IL[L] * Ylm(L, M, ax) * t2 * ap[L2M * num_lm * num_lm + l2mi * num_lm + l2mj];
                    }
                }
            }
            // Shift atom to the center instead of corner.
            for(int ix = 0; ix < dimx; ix++)
            {
                for(int iy = 0; iy < dimy; iy++)
                {
                    for(int iz = 0; iz < dimz; iz++)
                    {
                        int idx = ix * dimy * dimz + iy * dimz + iz;
                        std::complex<double> phaseshift =std::pow(phase, ix + ixstart + iy + iystart + iz + izstart);
                        betaptr[idx] *= phaseshift/vol;
                        betaptr_r[0][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);
                        betaptr_r[1][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);
                        betaptr_r[2][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);

                        // sqrt(3/4pi) is the normlaized constant for hamonic x, y, z and we only want non-normalized x, y, z
                    }
                }
            }

        }

        
        int amode = S_IREAD | S_IWRITE;
        std::string filename = "PROJECTORS/forward_beta_species" + std::to_string(isp) + "_kpt" + std::to_string(kindex);
        int fhand = open(filename.c_str(), O_CREAT | O_TRUNC | O_RDWR, amode);
        size_t count = sizeof(std::complex<double>) * sp->nh * pbasis;
        write(fhand, forward_beta, count);
        write(fhand, forward_beta_r[0], count);
        write(fhand, forward_beta_r[1], count);
        write(fhand, forward_beta_r[2], count);
        close(fhand);


    }  // end for

    delete [] IL;
    delete [] forward_beta;
    delete [] forward_beta_r[0];
    delete [] forward_beta_r[1];
    delete [] forward_beta_r[2];

} /* end InitDelocalizedWeight */

