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


// Used for localizing projectors in real space
template void Kpoint<double>::GetLocalizedWeight(void);
template void Kpoint<std::complex<double>>::GetLocalizedWeight(void);

template <class KpointType> void Kpoint<KpointType>::GetLocalizedWeight (void)
{

    int max_size;
    KpointType *Nlweight;
    std::complex<double> I_t(0.0, 1.0);

    int num_nonloc_ions = BetaProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = BetaProjector->get_nonloc_ions_list();


    SPECIES *sp;
    ION *iptr;

    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;
    std::complex<double> *beptr, *gbptr;
    std::complex<double> *in, *out;

    int P0_BASIS = get_P0_BASIS();

    /*maximum of nldim^3 for any species */
    max_size = ct.max_nldim * ct.max_nldim * ct.max_nldim;


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * 2 * max_size);
    std::complex<double> *phase_fftw = new std::complex<double>[max_size];

    if (beptr == NULL)
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    gbptr = beptr + max_size;


    Projector<KpointType> *P = BetaProjector;

    /* Loop over ions */
    for (int ion1 = 0; ion1 < num_nonloc_ions; ion1++)
    {

        int ion = nonloc_ions_list[ion1];
        /* Generate ion pointer */
        iptr = &Atoms[ion];

        /* Get species type */
        sp = &Species[iptr->species];

        Nlweight = &nl_weight[ion1 * ct.max_nl * P0_BASIS];

        int nlxdim = P->get_nldim(iptr->species);
        int nlydim = P->get_nldim(iptr->species);
        int nlzdim = P->get_nldim(iptr->species);

        in = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * nlxdim * nlydim * nlzdim);
        out = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * nlxdim * nlydim * nlzdim);


        /*Number of grid points on which fourier transform is done (in the coarse grid) */
        int coarse_size = nlxdim * nlydim * nlzdim;

        /*Calculate the phase factor */
        FindPhase(sp, nlxdim, nlydim, nlzdim, P->nlcrds[ion].data(), phase_fftw);

        /*Temporary pointer to the already calculated forward transform */
        fptr = (std::complex<double> *)&sp->forward_beta[kidx*sp->num_projectors * coarse_size];


        /* Loop over radial projectors */
        for (int ip = 0; ip < sp->num_projectors; ip++)
        {


            /*Apply the phase factor */
            for (int idx = 0; idx < coarse_size; idx++)
            {
                gbptr[idx] =fptr[idx] * std::conj(phase_fftw[idx]);
            }


            /*Do the backwards transform */
            sp->prj_pwave->FftInverse(gbptr, beptr);

            /*This takes and stores the part of beta that is useful for this PE */
            AssignWeight (this, sp, ion, reinterpret_cast<fftw_complex*>(beptr), Nlweight);


            /*Advance the temp pointers */
            fptr += coarse_size;
            Nlweight += P0_BASIS;

        }                   /*end for(ip = 0;ip < sp->num_projectors;ip++) */


        fftw_free(out);
        fftw_free(in);


    }                           /* end for */


    fftw_free (beptr);
    delete [] phase_fftw;


}                               /* end GetWeight */
