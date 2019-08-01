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

template void LcaoGetPsi(Kpoint<double> *);
template void LcaoGetPsi(Kpoint<std::complex<double> > *);

template <typename KpointType>
void LcaoGetPsi (Kpoint<KpointType> *kptr)
{

    State<KpointType> *states = kptr->Kstates;
    double *kvec = kptr->kp.kvec;

    long idum;

    int PX0_GRID = get_PX0_GRID();
    int PY0_GRID = get_PY0_GRID();
    int PZ0_GRID = get_PZ0_GRID();
    int PX_OFFSET = get_PX_OFFSET();
    int PY_OFFSET = get_PY_OFFSET();
    int PZ_OFFSET = get_PZ_OFFSET();
    int P0_BASIS = get_P0_BASIS();

    for (int st = 0; st < ct.num_states; st++)
    {
        for (int idx = 0; idx < P0_BASIS; idx++)
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

    if(state_count <= ct.num_states)
    {
        double coeff = 1.0;
        int st = 0;
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            /* Generate ion pointer */
            ION &Atom = Atoms[ion];

            /* Get species type */
            SPECIES &AtomType = Species[Atom.species];

            if(ct.atomic_orbital_type == DELOCALIZED)
            {
                // Delocalized orbitals
                kptr->get_ion_orbitals(&Atom, states[st].psi);
                st += AtomType.num_orbitals;
            }
            else
            {
                /*Loop over atomic wavefunctions for given ion*/
                for (int ip = 0; ip < AtomType.num_atomic_waves; ip++)
                {
                    int l = AtomType.atomic_wave_l[ip];
                    if(AtomType.atomic_wave_oc[ip] > 0.0)
                    {
                        /*Loop over all m values for given l and get wavefunctions */
                        for (int m=0; m < 2*l+1; m++)
                        {
                            LcaoGetAwave(states[st].psi, &Atom, ip, l, m, coeff, kvec);
                            st++;
                        }
                    }
                }
            }
        }
    }
    else
    {
        KpointType *npsi = new KpointType[P0_BASIS * state_count];      // The array for orbital storage is 4x run_state which should be big enough
        long *aidum = new long[state_count];
        for(int st = 0;st < state_count;st++)
        {
            aidum[st] = idum = st + 3314;
        }

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
                kptr->get_ion_orbitals(&Atom, &npsi[wave_idx * P0_BASIS]);
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
                            for(int idx = 0;idx < P0_BASIS;idx++)  npsi[wave_idx * P0_BASIS + idx] = 0.0;
                            LcaoGetAwave(&npsi[wave_idx * P0_BASIS], &Atom, ip, l, m, coeff, kvec);
                            wave_idx++;
                        }

                    }
                }
            }
        }


        // Now generate a random mix
#if GPU_ENABLED
        KpointType *rmatrix = (KpointType *)GpuMallocManaged(state_count * ct.num_states * sizeof(KpointType));
#else
        KpointType *rmatrix = new KpointType[state_count * ct.num_states];
#endif

        for(int st = 0;st < state_count;st++) {
            for(int idx = 0;idx < ct.num_states;idx++) {
                rmatrix[idx*state_count + st] = rand0(&aidum[idx]);
            }
        }

        char *trans_n = "n";
        KpointType alpha(1.0);
        KpointType beta(0.0);

    
        RmgGemm(trans_n, trans_n, P0_BASIS, ct.num_states, state_count, alpha,
            npsi, P0_BASIS, rmatrix, state_count, beta, states[0].psi, P0_BASIS);

#if GPU_ENABLED
        GpuFreeManaged(rmatrix);
#else
        delete [] rmatrix;
#endif


        delete [] npsi;
        delete [] aidum;

    }

    /*Initialize any additional states to random start*/
    if ( ct.num_states > state_count)
    {
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


        for (int st = state_count; st < ct.num_states; st++)
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
