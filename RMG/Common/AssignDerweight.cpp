/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "FiniteDiff.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

template void AssignDerweight<double> (Kpoint<double> *, SPECIES *, int, fftw_complex *, double *);
template void AssignDerweight<std::complex<double> >(Kpoint<std::complex<double>> *, SPECIES *, int, fftw_complex *, std::complex<double> *);


template <typename KpointType>
void AssignDerweight (Kpoint<KpointType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, KpointType *Nlweight)
{

    Lattice *L = kptr->L;
    TradeImages *T = kptr->T;


    int pbasis = kptr->pbasis;
    int nldim = sp->nldim;
    KpointType ZERO_t(0.0);

    std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight;

    double *Nlweight_R = (double *)Nlweight;

    weight_shift_center(sp, beptr);

    for(int idx = 0; idx < pbasis; idx++) Nlweight[idx] = ZERO_t;
    if(pct.idxptrlen[ion] == 0) return;

    std::complex<double> *nbeptr = (std::complex<double> *)beptr;


    double *pR = pct.phaseptr[ion];
    pR += 2 * kptr->kidx * nldim * nldim * nldim;
    double *pI = pR + nldim * nldim * nldim;

    int *pidx = pct.nlindex[ion];
    int *dvec = pct.idxflag[ion];
    int idx = 0;
    int docount = 0;


    idx = 0;
    docount = 0;
    for (int ix = 0; ix < sp->nldim; ix++)
    {

        for (int iy = 0; iy < sp->nldim; iy++)
        {

            for (int iz = 0; iz < sp->nldim; iz++)
            {

                if (dvec[idx])
                {
                    int idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                    if(ct.is_gamma) {
                        Nlweight_R[pidx[docount]] = std::real(nbeptr[idx1]);
                    }
                    else {
                        Nlweight_C[pidx[docount]] = std::complex<double>(std::real(nbeptr[idx1]) * pR[idx1], -std::real(nbeptr[idx1]) * pI[idx1]);
                    }
                    docount++;
                }

                idx++;
            }
        }
    }
    if (docount != pct.idxptrlen[ion])
    {
        rmg_printf ("docount = %d != %d = pct.idxptrlen[ion = %d]\n", docount, pct.idxptrlen[ion], ion);
        rmg_error_handler (__FILE__, __LINE__, "wrong numbers of projectors");
    }


}
