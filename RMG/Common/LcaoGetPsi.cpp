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

    if(ct.spinorbit && state_count > nstates)
    {
        rmg_printf("state_count %d != nstates %d", state_count, nstates);
        rmg_error_handler(__FILE__,__LINE__," state_count != nstates Terminating.");

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

    if(ct.spinorbit || ct.noncoll)
    {
        wave_idx = 0;
        int state_idx = 0;
        int tot_LM = (ct.max_l +1) *(ct.max_l +1);
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            ION &Atom = Atoms[ion];
            SPECIES &AtomType = Species[Atom.species];
            if(AtomType.is_spinorb)
            {
                std::complex<double> *Umm = AtomType.Umm;
                for (int ip = 0; ip < AtomType.num_atomic_waves; ip++)
                {
                    int li = AtomType.atomic_wave_l[ip];
                    double ji = AtomType.atomic_wave_j[ip];

                    if( std::abs(ji -li - 0.5) < 1.0e-5)
                    {
                        for(int m = -li -1; m <= li; m++)
                        {
                            double alpha_up = std::sqrt( (li + m + 1.0)/(2*li + 1.0));
                            double alpha_dn = std::sqrt( (li - m )/(2*li + 1.0));
                            int lmm_up = li * li + li-m;
                            int lmm_dn = li * li + li - (m+1);
                            std::complex<double> *psi_C = (std::complex<double> *)states[state_idx].psi;

                            for(int mp = 0; mp < 2*li+1; mp++)
                            {

                                int lmp = li * li + mp;

                                for(int idx = 0; idx < pbasis; idx++)
                                {
                                    if(m  >=-li) 
                                        psi_C[idx] += alpha_up * Umm[lmm_up * tot_LM + lmp] *  npsi[(wave_idx+mp) * pbasis + idx];
                                    if(m  < li) 
                                        psi_C[idx + pbasis] += alpha_dn * Umm[lmm_dn * tot_LM + lmp] *  npsi[(wave_idx+mp) * pbasis + idx];
                                }
                            }
                            state_idx++;
                        }
                    }
                    //  case: j = l - 1/2, mj =[-j,j], m = mj+1/2 = [-l+1, l]
                    else if( std::abs(ji -li + 0.5) < 1.0e-5)
                    {
                        for(int m = -li+1; m <= li; m++)
                        {
                            double alpha_up = std::sqrt( (li - m + 1.0)/(2*li + 1.0));
                            double alpha_dn = -std::sqrt( (li + m )/(2*li + 1.0));
                            int lmm_up = li * li + li-(m-1);
                            int lmm_dn = li * li + li-m;
                            std::complex<double> *psi_C = (std::complex<double> *)states[state_idx].psi;

                            for(int mp = 0; mp < 2*li+1; mp++)
                            {

                                int lmp = li * li + mp;

                                for(int idx = 0; idx < pbasis; idx++)
                                {
                                    psi_C[idx] += alpha_up * Umm[lmm_up * tot_LM + lmp] *  npsi[(wave_idx+mp) * pbasis + idx];
                                    psi_C[idx + pbasis] += alpha_dn * Umm[lmm_dn * tot_LM + lmp] *  npsi[(wave_idx+mp) * pbasis + idx];
                                }
                            }
                            state_idx++;

                        }
                    }

                    wave_idx += 2 * li + 1;
                }
            }
            else
            {
                for (int ip = 0; ip < AtomType.num_atomic_waves; ip++)
                {
                    int l = AtomType.atomic_wave_l[ip];
                    for (int m=0; m < 2*l+1; m++)
                    {
                        for(int idx = 0; idx < pbasis; idx++)
                        {
                            states[state_idx   ].psi[idx] += npsi[wave_idx * pbasis + idx];
                            states[state_idx +1].psi[idx+pbasis] += npsi[wave_idx * pbasis + idx];
                        }

                        state_idx += 2;
                        wave_idx ++;
                    }
                }
            }


        }
    }
    else if(state_count <= nstates)
    {
        for(int st = 0;st < wave_idx;st++)
        {
            for(int idx = 0; idx < pbasis; idx++)
                states[st].psi[idx] = npsi[st * pbasis + idx];
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
    if ( nstates > state_count)
    {
        int ix, iy, iz;
        int xoff, yoff, zoff;
        State<KpointType> *state_p;

        int NX_GRID = get_NX_GRID();
        int NY_GRID = get_NY_GRID();
        int NZ_GRID = get_NZ_GRID();
        double *xrand = new double[2 * NX_GRID];
        double *yrand = new double[2 * NY_GRID];
        double *zrand = new double[2 * NZ_GRID];

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
            for (int idx = 0; idx < 2*NX_GRID; idx++)
                xrand[idx] = rand0 (&idum) - 0.5;
            for (int idx = 0; idx < 2*NY_GRID; idx++)
                yrand[idx] = rand0 (&idum) - 0.5;
            for (int idx = 0; idx < 2*NZ_GRID; idx++)
                zrand[idx] = rand0 (&idum) - 0.5;

            state_p = &states[st];

            int idx = 0;
            for (int ix = 0; ix < PX0_GRID; ix++)
            {

                for (int iy = 0; iy < PY0_GRID; iy++)
                {

                    for (int iz = 0; iz < PZ0_GRID; iz++)
                    {

                        for(int is = 0; is < ct.noncoll_factor; is++)
                        {
                            state_p->psi[idx + is * pbasis] = xrand[is * NX_GRID + xoff + ix] * 
                                yrand[is * NY_GRID + yoff + iy] * zrand[is * NZ_GRID + zoff + iz];
                            state_p->psi[idx + is * pbasis] = state_p->psi[idx + is * pbasis] * state_p->psi[idx + is * pbasis];
                        }
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

