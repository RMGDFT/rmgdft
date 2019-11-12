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

#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "transition.h"

template void Kpoint<double>::LcaoGetPsi(void);
template void Kpoint<std::complex<double>>::LcaoGetPsi(void);

template <class KpointType> void Kpoint<KpointType>::LcaoGetPsi (void)
{

    State<KpointType> *states = Kstates;
    double *kvec = kp.kvec;

    long idum;

    int PX0_GRID = get_PX0_GRID();
    int PY0_GRID = get_PY0_GRID();
    int PZ0_GRID = get_PZ0_GRID();
    int PX_OFFSET = get_PX_OFFSET();
    int PY_OFFSET = get_PY_OFFSET();
    int PZ_OFFSET = get_PZ_OFFSET();

    for (int st = 0; st < nstates; st++)
    {
        for (int idx = 0; idx < pbasis * ct.noncoll_factor; idx++)
        {
            states[st].psi[idx] = 0.0;
        }
    }

    idum = 1314; 
    rand0 (&idum);

    /* Loop over ions */
    int state_count = 0;

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        /* Get species type */
        SPECIES &AtomType = Species[Atoms[ion].species];

        /*Make sure that the wavefunctions have been read*/
        if (!AtomType.num_atomic_waves) {
            rmg_printf("No initial wavefunctions for ion %lu, most likely the PP file does not have them", ion);
            rmg_error_handler(__FILE__,__LINE__,"Terminating.");
        }
    }

    state_count = CountAtomicOrbitals();

    KpointType *npsi = (KpointType *)GpuMallocManaged(state_count * pbasis * sizeof(KpointType));

    double coeff = 1.0;
    int wave_idx = 0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        /* Generate ion pointer */
        ION &Atom = Atoms[ion];

        /* Get species type */
        SPECIES &AtomType = Species[Atom.species];

        if(ct.atomic_orbital_type == DELOCALIZED)
        {
            // Delocalized orbitals
            get_ion_orbitals(&Atom, &npsi[wave_idx * pbasis]);
            wave_idx += AtomType.num_orbitals;
        }
        else
        {
            /*Loop over atomic wavefunctions for given ion*/
            for (int ip = 0; ip < AtomType.num_atomic_waves; ip++)
            {
                int l = AtomType.atomic_wave_l[ip];
                if(AtomType.atomic_wave_oc[ip] > 0.0) {

                    /*Loop over all m values for given l and get wavefunctions */
                    for (int m=0; m < 2*l+1; m++)
                    {
                        for(int idx = 0;idx < pbasis;idx++)  npsi[wave_idx * pbasis + idx] = 0.0;
                        LcaoGetAwave(&npsi[wave_idx * pbasis], &Atom, ip, l, m, coeff, kvec);
                        wave_idx++;
                    }

                }
            }
        }
    }

    // in the case of noncollinear, the first state_count wavefunctions will have (atomic_psi, 0) and the second state_count wavefuntions will
    // have (0, atomic_psi). 
    std::cout<< "lcao " << ct.noncoll_factor << std::endl;
    if(state_count * ct.noncoll_factor <= nstates)
    {
        for(int st = 0;st < state_count;st++)
        {
            for(int idx = 0; idx < pbasis; idx++)
                states[st].psi[idx] = npsi[st * pbasis + idx];

            if(ct.noncoll)
                for(int idx = 0; idx < pbasis; idx++)
                    states[st+state_count].psi[idx+pbasis] = npsi[st * pbasis + idx];

        }
    }
    else
    {
        long *aidum = new long[state_count];
        for(int st = 0;st < state_count;st++)
        {
            aidum[st] = idum = st + 3314;
        }

        // Now generate a random mix
        KpointType *rmatrix = (KpointType *)GpuMallocManaged(state_count * nstates * sizeof(KpointType));

        for(int st = 0;st < state_count;st++) {
            for(int idx = 0;idx < nstates;idx++) {
                rmatrix[idx*state_count + st] = rand0(&aidum[idx]);
            }
        }

        char *trans_n = "n";
        KpointType alpha(1.0);
        KpointType beta(0.0);


        int lda = pbasis * ct.noncoll_factor;
        RmgGemm(trans_n, trans_n, pbasis, nstates, state_count, alpha,
                npsi, pbasis, rmatrix, state_count, beta, states[0].psi, lda);

        if(ct.noncoll)
            RmgGemm(trans_n, trans_n, pbasis, nstates, state_count, alpha,
                    npsi, pbasis, rmatrix, state_count, beta, states[state_count].psi+pbasis, lda);


        GpuFreeManaged(rmatrix);
        delete [] aidum;
    }
    GpuFreeManaged(npsi);


    /*Initialize any additional states to random start*/
    if ( nstates > state_count * ct.noncoll_factor)
    {
        if(ct.noncoll) 
        {
            rmg_printf("in the case of noncollinear, nstates cannot be bigger than number of atomic wavefunctions\n");
            rmg_error_handler(__FILE__,__LINE__,"Terminating.");
        }
        int ix, iy, iz;
        int xoff, yoff, zoff;
        State<KpointType> *state_p;

        double *xrand = new double[2 * get_NX_GRID()];
        double *yrand = new double[2 * get_NY_GRID()];
        double *zrand = new double[2 * get_NZ_GRID()];

        pe2xyz (pct.gridpe, &ix, &iy, &iz);
        xoff = PX_OFFSET;
        yoff = PY_OFFSET;
        zoff = PZ_OFFSET;

        /* Initialize the random number generator */
        idum = 3356;
        rand0 (&idum);


        for (int st = state_count; st < nstates; st++)
        {

            /* Generate x, y, z random number sequences */
            for (int idx = 0; idx < get_NX_GRID(); idx++)
                xrand[idx] = rand0 (&idum) - 0.5;
            for (int idx = 0; idx < get_NY_GRID(); idx++)
                yrand[idx] = rand0 (&idum) - 0.5;
            for (int idx = 0; idx < get_NZ_GRID(); idx++)
                zrand[idx] = rand0 (&idum) - 0.5;

            state_p = &states[st];

            int idx = 0;
            for (int ix = 0; ix < PX0_GRID; ix++)
            {

                for (int iy = 0; iy < PY0_GRID; iy++)
                {

                    for (int iz = 0; iz < PZ0_GRID; iz++)
                    {

                        state_p->psi[idx] = xrand[xoff + ix] * yrand[yoff + iy] * zrand[zoff + iz];
                        state_p->psi[idx] = state_p->psi[idx] * state_p->psi[idx];
                        idx++;

                    }               /* end for */
                }                   /* end for */
            }                       /* end for */


        }                           /* end for */

        delete [] zrand;
        delete [] yrand;
        delete [] xrand;

    }


}
/******/
