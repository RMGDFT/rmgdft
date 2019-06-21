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


// Used to generate orbital projectors that span the full space
template void GetDelocalizedOrbital<double> (Kpoint<double> **Kptr);
template void GetDelocalizedOrbital<std::complex<double> >(Kpoint<std::complex<double>> **Kptr);
template <typename KpointType>
void GetDelocalizedOrbital (Kpoint<KpointType> **Kptr)
{

    KpointType *weight;
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
    KpointType *tem_array = new KpointType[pbasis];

    for(int kpt =0; kpt < ct.num_kpts_pe;kpt++) {

        double *kvec = Kptr[kpt]->kvec;

        Projector<KpointType> *P = Kptr[kpt]->OrbitalProjector;
        size_t stride = P->get_pstride();

        /* Loop over ions */
        for (int ion = 0; ion < ct.num_ions; ion++)
        {
            size_t offset = (size_t)ion * stride * (size_t)pbasis;
            weight = &Kptr[kpt]->orbital_weight[offset];

            /* Generate ion pointer */
            ION *iptr = &ct.ions[ion];

            /* Get species type */
            SPECIES *sp = &ct.sp[iptr->species];

            int nlxdim = get_NX_GRID();
            int nlydim = get_NY_GRID();
            int nlzdim = get_NZ_GRID();

            /*Calculate the phase factor for delocalized case */
            FindPhaseKpoint (kvec, nlxdim, nlydim, nlzdim, P->nlcrds[ion].data(), fftw_phase, false);

            /*Temporary pointer to the already calculated forward transform */
            fptr = (std::complex<double> *)&sp->forward_orbital[kpt * sp->num_orbitals * pbasis];


            /* Loop over radial projectors */
            for (int ip = 0; ip < sp->num_ldaU_orbitals; ip++)
            {


                /*Apply the phase factor */
                for (int idx = 0; idx < pbasis; idx++)
                {
                    gbptr[idx] =  fptr[idx] * std::conj(fftw_phase[idx]);
                }


                /*Do the backwards transform */
                PfftInverse(gbptr, beptr, *coarse_pwaves);

                std::complex<double> *weight_C = (std::complex<double> *)weight;

                double *weight_R = (double *)weight;

                for(int idx = 0; idx < pbasis; idx++) weight[idx] = ZERO_t;

                std::complex<double> *nbeptr = (std::complex<double> *)beptr;

                std::complex<double> *tem_array_C = (std::complex<double> *)tem_array;

                if(ct.is_gamma)
                {
                    for(int ix = 0; ix <pbasis; ix++) tem_array[ix] = std::real(nbeptr[ix]);
                }
                else
                {
                    for (int idx = 0; idx < pbasis; idx++) tem_array_C[idx] = nbeptr[idx];
                }

                for (int idx = 0; idx < pbasis; idx++)
                {
                    if(ct.is_gamma) {
                        weight_R[idx] = std::real(nbeptr[idx]);

                    }
                    else {
                        weight_C[idx] = nbeptr[idx];
                    }
                }


                /*Advance the temp pointers */
                fptr += pbasis;
                weight += pbasis;

            }                   /*end for(ip = 0;ip < sp->num_ldaU_orbitals;ip++) */


        }                           /* end for */

    } // end for(kpt)

    delete [] tem_array;
    delete [] fftw_phase;
    fftw_free (gbptr);
    fftw_free (beptr);


}                               /* end GetWeight */
