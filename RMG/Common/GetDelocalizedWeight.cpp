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
template void Kpoint<double>::GetDelocalizedWeight(void);
template void Kpoint<std::complex<double>>::GetDelocalizedWeight(void);

template <class KpointType> void Kpoint<KpointType>::GetDelocalizedWeight (void)
{

    KpointType *Nlweight;
    KpointType ZERO_t(0.0);
    std::complex<double> I_t(0.0, 1.0);

    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);

    if ((beptr == NULL) || (gbptr == NULL))
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];

    double *kvec = kp.kvec;

    Projector<KpointType> *P = BetaProjector;
    size_t stride = P->get_pstride();

    /* Loop over ions */
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        size_t offset = (size_t)ion * stride * (size_t)pbasis;

        /* Get species type */
        SPECIES &AtomType = Species[Atoms[ion].species];


        int nlxdim = get_NX_GRID();
        int nlydim = get_NY_GRID();
        int nlzdim = get_NZ_GRID();

        /*Calculate the phase factor for delocalized case */

        FindPhaseKpoint (kvec, nlxdim, nlydim, nlzdim, P->nlcrds[ion].data(), fftw_phase, false);

        /* Loop over radial projectors */
        for (int ip = 0; ip < AtomType.num_projectors; ip++)
        {
            Nlweight = &nl_weight[offset + ip * pbasis];

            /*Temporary pointer to the already calculated forward transform */
            fptr = (std::complex<double> *)AtomType.forward_beta[kidx * AtomType.num_projectors * pbasis + ip*pbasis];

            /*Apply the phase factor */
            for (int idx = 0; idx < pbasis; idx++) gbptr[idx] =  fptr[idx] * std::conj(fftw_phase[idx]);

            /*Do the backwards transform */
            coarse_pwaves->FftInverse(gbptr, beptr);

            std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight;
            double *Nlweight_R = (double *)Nlweight;
            for(int idx = 0; idx < pbasis; idx++) Nlweight[idx] = ZERO_t;

            std::complex<double> *nbeptr = (std::complex<double> *)beptr;

            // Apply B operator then map weights back
            if(ct.is_gamma) {
                for (int idx = 0; idx < pbasis; idx++) Nlweight_R[idx] = std::real(nbeptr[idx]);
            }
            else {
                for (int idx = 0; idx < pbasis; idx++) Nlweight_C[idx] = nbeptr[idx];
            }
        }

        // for stress calculation, calculate x*beta, y * beta, z*beta
        if(!ct.stress) continue;
        for(int ixyz = 0; ixyz < 3; ixyz++)
        {
            for (int ip = 0; ip < AtomType.num_projectors; ip++)
            {
                Nlweight = &nl_weight[nl_weight_size * (ixyz+1) + offset + ip * pbasis];

                /*Temporary pointer to the already calculated forward transform */
                fptr = (std::complex<double> *)AtomType.forward_beta_r[ixyz][kidx * AtomType.num_projectors * pbasis + ip*pbasis];

                /*Apply the phase factor */
                for (int idx = 0; idx < pbasis; idx++) gbptr[idx] =  fptr[idx] * std::conj(fftw_phase[idx]);

                /*Do the backwards transform */
                coarse_pwaves->FftInverse(gbptr, beptr);

                std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight;
                double *Nlweight_R = (double *)Nlweight;
                for(int idx = 0; idx < pbasis; idx++) Nlweight[idx] = ZERO_t;

                std::complex<double> *nbeptr = (std::complex<double> *)beptr;

                // Apply B operator then map weights back
                if(ct.is_gamma) {
                    for (int idx = 0; idx < pbasis; idx++) Nlweight_R[idx] = std::real(nbeptr[idx]);
                }
                else {
                    for (int idx = 0; idx < pbasis; idx++) Nlweight_C[idx] = nbeptr[idx];
                }

            } 
        }

    }                           /* end for */



    delete [] fftw_phase;
    fftw_free (gbptr);
    fftw_free (beptr);

#if HIP_ENABLED
    size_t stress_factor = 1;
    if(ct.stress) stress_factor = 4;
    hipMemcpy(nl_weight_gpu, nl_weight, stress_factor*nl_weight_size*sizeof(KpointType), hipMemcpyHostToDevice);
#elif CUDA_ENABLED
    size_t stress_factor = 1;
    if(ct.stress) stress_factor = 4;
    cudaMemcpy(nl_weight_gpu, nl_weight, stress_factor*nl_weight_size*sizeof(KpointType), cudaMemcpyHostToDevice);
#endif

}                               /* end GetWeight */
