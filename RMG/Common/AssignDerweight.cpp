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


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"


template void AssignDerweight<double> (Kpoint<double> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, double * rtptr);
template void AssignDerweight<std::complex<double> > (Kpoint<std::complex<double>> *kptr, SPECIES * sp, int ion, fftw_complex * beptr,
std::complex<double> *rtptr);



template <typename OrbitalType>
void AssignDerweight (Kpoint<OrbitalType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, OrbitalType *rtptr)
        
{

    int idx, ix, iy, iz, *dvec;
    int idx1, docount;

    weight_shift_center(sp, beptr);
    int *pidx = pct.nlindex[ion];

    double *pR = pct.phaseptr[ion];
    pR += 2 * kptr->kidx * sp->nldim * sp->nldim *sp->nldim;
    double *pI = pR + sp->nldim * sp->nldim * sp->nldim;

    double *rtptr_R = (double *)rtptr;
    std::complex<double> *rtptr_C = (std::complex<double> *)rtptr;

    std::complex<double> *nbeptr = (std::complex<double> *)beptr;


    dvec = pct.idxflag[ion];
    idx = docount = 0;
    for (ix = 0; ix < sp->nldim; ix++)
    {

        for (iy = 0; iy < sp->nldim; iy++)
        {

            for (iz = 0; iz < sp->nldim; iz++)
            {

                if (dvec[idx])
                {
                    idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;

                    if (std::imag(nbeptr[idx1]) > 1.0e-8)
                    {
                        printf ("beptr[%d].im=%e\n", idx1, std::imag(nbeptr[idx1]));
                        rmg_error_handler (__FILE__, __LINE__, "something wrong with the fourier transformation");
                    }

                    if(ct.is_gamma)
                        rtptr_R[pidx[docount]] = std::real(nbeptr[idx1]);
                    else
                        rtptr_C[pidx[docount]] = std::real(nbeptr[idx1]) * std::complex<double>(pR[idx1], -pI[idx1]);
                    docount++;
                }

                idx++;
            }
        }
    }
    if (docount != pct.idxptrlen[ion])
    {
        printf ("docount = %d != %d = pct.idxptrlen[ion = %d]\n", docount, pct.idxptrlen[ion], ion);
        rmg_error_handler (__FILE__, __LINE__, "wrong numbers of projectors");
    }
}
