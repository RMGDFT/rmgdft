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
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgParallelFft.h"


// Used to generate projectors that span the full space
template void GetWeight<double> (Kpoint<double> **Kptr);
template void GetWeight<std::complex<double> >(Kpoint<std::complex<double>> **Kptr);
template <typename KpointType>
void GetWeight (Kpoint<KpointType> **Kptr)
{

    int max_size;
    double *rtptr;
    KpointType *Bweight, *Nlweight;
    KpointType ZERO_t(0.0);
    std::complex<double> I_t(0.0, 1.0);
    int pbasis = Kptr[0]->pbasis;

    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;
    std::complex<double> *beptr, *gbptr;

    /*maximum of nldim^3 for any species */
    max_size = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * 2 * max_size);

    if (beptr == NULL)
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    gbptr = beptr + max_size;

    // Release memory allocated for fftw_phase_sin and fftw_phase_cos prior to 
    // reallocation (if needed) in find_phase
    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];

    for(int idx = 0; idx < pct.num_tot_proj * pbasis; idx++)
    {
        pct.weight[idx] = 0.0;
        pct.Bweight[idx] = 0.0;
    }

    for(int kpt =0; kpt < ct.num_kpts_pe;kpt++) {

        KpointType *tem_array = new KpointType[pbasis];
        KpointType *Btem_array = new KpointType[pbasis];

        /* Loop over ions */
        for (int ion = 0; ion < ct.num_ions; ion++)
        {
            rtptr = &pct.weight[ion * ct.max_nl * pbasis];
            Bweight = &Kptr[kpt]->nl_Bweight[ion * ct.max_nl * pbasis];
            Nlweight = &Kptr[kpt]->nl_weight[ion * ct.max_nl * pbasis];

            /* Generate ion pointer */
            ION *iptr = &ct.ions[ion];

            /* Get species type */
            SPECIES *sp = &ct.sp[iptr->species];

            int nlxdim = get_NX_GRID();
            int nlydim = get_NY_GRID();
            int nlzdim = get_NZ_GRID();

            /*Calculate the phase factor */
            FindPhase (nlxdim, nlydim, nlzdim, iptr->nlcrds, fftw_phase);

            /*Temporary pointer to the already calculated forward transform */
            fptr = (std::complex<double> *)&sp->forward_beta[kpt * sp->num_projectors * pbasis];


            /* Loop over radial projectors */
            for (int ip = 0; ip < sp->num_projectors; ip++)
            {


                /*Apply the phase factor */
                for (int idx = 0; idx < pbasis; idx++)
                {
                    gbptr[idx] =  fptr[idx] * std::conj(fftw_phase[idx]);
                }


                /*Do the backwards transform */
                PfftInverse(gbptr, beptr, *coarse_pwaves);

                std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight;
                std::complex<double> *Bweight_C = (std::complex<double> *)Bweight;

                double *Nlweight_R = (double *)Nlweight;
                double *Bweight_R = (double *)Bweight;

                for(int idx = 0; idx < pbasis; idx++) rtptr[idx] = 0.0;
                for(int idx = 0; idx < pbasis; idx++) Bweight[idx] = ZERO_t;
                for(int idx = 0; idx < pbasis; idx++) Nlweight[idx] = ZERO_t;

                std::complex<double> *nbeptr = (std::complex<double> *)beptr;

                std::complex<double> *tem_array_C = (std::complex<double> *)tem_array;
                std::complex<double> *Btem_array_C = (std::complex<double> *)Btem_array;
                double *Btem_array_R = (double *)Btem_array;


                for(int ix = 0; ix <pbasis; ix++) {
                    tem_array[ix] = std::real(nbeptr[ix]);
                }


                // Apply phase factor for non-gamma
// FIX: have to fix this up for non-gamma
                if(!ct.is_gamma) {
                    for (int idx = 0; idx < pbasis; idx++)
                        tem_array_C[idx] = nbeptr[idx];
                }

                // Apply B operator then map weights back
//                AppCirDriverBeta (Kptr[kpt]->L, Kptr[kpt]->T, tem_array, Btem_array, nlxdim, nlydim, nlzdim, ct.kohn_sham_fd_order);
// FIX: Only works for central operator right now. Have to fix AppCirDriverBeta for Mehrstellen
for(int idx = 0;idx < pbasis;idx++)Btem_array[idx] = tem_array[idx];

                for (int idx = 0; idx < pbasis; idx++)
                {
                    if(ct.is_gamma) {
                        Nlweight_R[idx] = std::real(nbeptr[idx]);
                        Bweight_R[idx] = Btem_array_R[idx];

                    }
                    else {
                        rtptr[idx] = std::real(nbeptr[idx]);
                        Nlweight_C[idx] = nbeptr[idx];
                        Bweight_C[idx] = Btem_array_C[idx];
                    }

                }


                /*Advance the temp pointers */
                fptr += pbasis;
                rtptr += pbasis;
                Bweight += pbasis;
                Nlweight += pbasis;

            }                   /*end for(ip = 0;ip < sp->num_projectors;ip++) */


        }                           /* end for */

        delete [] Btem_array;
        delete [] tem_array;


    } // end for(kpt)

    delete [] fftw_phase;
    fftw_free (beptr);


}                               /* end get_weight */
