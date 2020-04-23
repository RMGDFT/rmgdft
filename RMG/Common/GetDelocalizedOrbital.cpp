/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
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


// Used to generate LDA+U orbital projectors that span the full space
template void Kpoint<double>::GetDelocalizedOrbital(void);
template void Kpoint<std::complex<double>>::GetDelocalizedOrbital(void);
template <class KpointType> void Kpoint<KpointType>::GetDelocalizedOrbital (void)
{

    KpointType *weight;
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

        Projector<KpointType> *P = OrbitalProjector;
        size_t stride = P->get_pstride();

        /* Loop over ions */
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            size_t offset = (size_t)ion * stride * (size_t)pbasis;
            weight = &orbital_weight[offset];
    
            /* Generate atom reference */
            ION &Atom = Atoms[ion];

            /* Get species type */
            SPECIES *sp = &Species[Atom.species];

            int nlxdim = get_NX_GRID();
            int nlydim = get_NY_GRID();
            int nlzdim = get_NZ_GRID();

            /*Calculate the phase factor for delocalized case */
            FindPhaseKpoint (kvec, nlxdim, nlydim, nlzdim, P->nlcrds[ion].data(), fftw_phase, false);

            /*Temporary pointer to the already calculated forward transform */
            fptr = (std::complex<double> *)&sp->forward_orbital[kidx * sp->num_orbitals * pbasis];


            /* Loop over radial projectors */
            for (int ip = 0; ip < sp->num_orbitals; ip++)
            {

                // This ranges over all orbitals including the m-dependence
                if(sp->awave_is_ldaU[ip])
                {

                    // Apply the phase factor.
                    for (int idx = 0; idx < pbasis; idx++) gbptr[idx] =  fptr[idx] * std::conj(fftw_phase[idx]);

                    /*Do the backwards transform */
                    coarse_pwaves->FftInverse(gbptr, beptr);

                    int lm = sp->awave_ldaU_lm[ip];
                    double awave_j = sp->awave_ldaU_j[ip];
                    double factor_j = 1.0;
                    if(sp->ldaU_l != 0 && std::abs(awave_j) > 1.0e-4)
                    {
                        factor_j = (2.0 * awave_j + 1.0) /(2.0 * sp->ldaU_l +1 )/2.0;
                    }

                    if(ct.is_gamma)
                    {
                        double *weight_R = (double *)&weight[lm * pbasis];
                        for (int idx = 0; idx < pbasis; idx++) weight_R[idx] += std::real(beptr[idx]) * factor_j;
                    }
                    else
                    {
                        std::complex<double> *weight_C = (std::complex<double> *)&weight[lm*pbasis];
                        for (int idx = 0; idx < pbasis; idx++) weight_C[idx] += beptr[idx] * factor_j;
                    }


                }

                /*Advance the temp pointers */
                fptr += pbasis;

            } /* end for(ip = 0;ip < sp->num_orbitals;ip++) */


        }                           /* end for */


    delete [] fftw_phase;
    fftw_free (gbptr);
    fftw_free (beptr);


}                               /* end GetWeight */
