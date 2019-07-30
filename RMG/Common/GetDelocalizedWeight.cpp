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


#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgParallelFft.h"


// Used to generate projectors that span the full space
template void GetDelocalizedWeight<double> (Kpoint<double> **Kptr);
template void GetDelocalizedWeight<std::complex<double> >(Kpoint<std::complex<double>> **Kptr);
template <typename KpointType>
void GetDelocalizedWeight (Kpoint<KpointType> **Kptr)
{

    KpointType *Nlweight;
    KpointType ZERO_t(0.0);
    std::complex<double> I_t(0.0, 1.0);
    int pbasis = Kptr[0]->pbasis;

    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);

    if ((beptr == NULL) || (gbptr == NULL))
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];

    for(int kpt =0; kpt < ct.num_kpts_pe;kpt++) {

        double *kvec = Kptr[kpt]->kvec;

        Projector<KpointType> *P = Kptr[kpt]->BetaProjector;
        size_t stride = P->get_pstride();

        /* Loop over ions */
        for (int ion = 0; ion < ct.num_ions; ion++)
        {
            size_t offset = (size_t)ion * stride * (size_t)pbasis;
            Nlweight = &Kptr[kpt]->nl_weight[offset];

            /* Generate ion pointer */
            ION *iptr = &Atoms[ion];

            /* Get species type */
            SPECIES *sp = &ct.sp[iptr->species];

            int nlxdim = get_NX_GRID();
            int nlydim = get_NY_GRID();
            int nlzdim = get_NZ_GRID();

            /*Calculate the phase factor for delocalized case */
            FindPhaseKpoint (kvec, nlxdim, nlydim, nlzdim, P->nlcrds[ion].data(), fftw_phase, false);

            /*Temporary pointer to the already calculated forward transform */
            fptr = (std::complex<double> *)&sp->forward_beta[kpt * sp->num_projectors * pbasis];


            /* Loop over radial projectors */
            for (int ip = 0; ip < sp->num_projectors; ip++)
            {

                /*Apply the phase factor */
                for (int idx = 0; idx < pbasis; idx++) gbptr[idx] =  fptr[idx] * std::conj(fftw_phase[idx]);

                /*Do the backwards transform */
                PfftInverse(gbptr, beptr, *coarse_pwaves);

                std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight;

                double *Nlweight_R = (double *)Nlweight;
                for(int idx = 0; idx < pbasis; idx++) Nlweight[idx] = ZERO_t;

                std::complex<double> *nbeptr = (std::complex<double> *)beptr;

// Apply B operator then map weights back
// FIX: Only works for central operator right now. Have to fix AppCirDriverBeta for Mehrstellen
                if(ct.is_gamma) {
                    for (int idx = 0; idx < pbasis; idx++) Nlweight_R[idx] = std::real(nbeptr[idx]);
                }
                else {
                    for (int idx = 0; idx < pbasis; idx++) Nlweight_C[idx] = nbeptr[idx];
                }

                /*Advance the temp pointers */
                fptr += pbasis;
                Nlweight += pbasis;

            }                   /*end for(ip = 0;ip < sp->num_projectors;ip++) */

        }                           /* end for */


    } // end for(kpt)

    delete [] fftw_phase;
    fftw_free (gbptr);
    fftw_free (beptr);


}                               /* end GetWeight */
