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

    // State array is always at least 2*ct.init_states so we can use the upper part for temp storage
    KpointType *npsi = Kstates[ct.init_states].psi;

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
                    if( !(AtomType.atomic_wave_oc[ip] > 0.0) ) continue;
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
                    if( !(AtomType.atomic_wave_oc[ip] > 0.0) ) continue;
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
        KpointType *rmatrix = (KpointType *)RmgMallocHost(state_count * nstates * sizeof(KpointType));

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


        RmgFreeHost(rmatrix);
        delete [] aidum;
    }


    /*Initialize any additional states to random start*/
    if ( nstates > state_count)
    {

        int NX_GRID = get_NX_GRID();
        int NY_GRID = get_NY_GRID();
        int NZ_GRID = get_NZ_GRID();
        BaseGrid *LG = new BaseGrid(Rmg_G->get_NX_GRID(1), Rmg_G->get_NY_GRID(1), Rmg_G->get_NZ_GRID(1), 1, 1, 1, 0, 1);
        int rank = Rmg_G->get_rank();
        MPI_Comm lcomm;
        MPI_Comm_split(Rmg_G->comm, rank+1, rank, &lcomm);
        LG->set_rank(0, lcomm);
        Pw *pwave = new Pw(*LG, Rmg_L, 1, false);
        double tpiba = 2.0*PI / Rmg_L.celldm[0];
        double crds[3], xtal[3];
        int half_size = pwave->pbasis;
        std::vector<int> index_sort(half_size);
        std::vector<int>gm(half_size);


        for(int ix = 0; ix<NX_GRID; ix++){
            for(int iy = 0; iy<NY_GRID; iy++){
                for(int iz = 0; iz<NZ_GRID; iz++){
                    int idx = ix * NZ_GRID * NY_GRID + iy * NZ_GRID + iz;
                    gm[idx] = 0;
                    if(iz >= NZ_GRID/2) gm[idx]  = 1;
                    if(iz == 0 && iy > NY_GRID/2) gm[idx]  = 1;
                    if(iz == 0 && iy == 0 && ix > NX_GRID/2) gm[idx]  = 1;

                }
            }
        }


        for(int idx = 0; idx < half_size; idx++) index_sort[idx] = idx;
        std::stable_sort(index_sort.begin(), index_sort.end(), [&](int i1, int i2) { return pwave->gmags[i1] < pwave->gmags[i2]; } );   

        int xoff, yoff, zoff;


        xoff = PX_OFFSET;
        yoff = PY_OFFSET;
        zoff = PZ_OFFSET;


        for (int st = state_count; st < nstates; st+=ct.noncoll_factor)
        {

            int idx = 0;
            int idx_pwave = index_sort[st - state_count/ct.noncoll_factor  ];
            for (int ix = 0; ix < PX0_GRID; ix++)
            {

                for (int iy = 0; iy < PY0_GRID; iy++)
                {

                    for (int iz = 0; iz < PZ0_GRID; iz++)
                    {
                        xtal[0] = (xoff+ix)/(double)NX_GRID;
                        xtal[1] = (yoff+iy)/(double)NY_GRID;
                        xtal[2] = (zoff+iz)/(double)NZ_GRID;
                        Rmg_L.to_cartesian(xtal, crds);
                        double theta = crds[0] * (tpiba*pwave->g[idx_pwave].a[0] ) +
                            crds[1] * (tpiba*pwave->g[idx_pwave].a[1] ) +
                            crds[2] * (tpiba*pwave->g[idx_pwave].a[2] );

                        if(gm[idx_pwave] == 1)
                        {
                            states[st].psi[idx] = sin(theta);
                        }
                        else
                        {
                            states[st].psi[idx] = cos(theta);
                        }
                        if(ct.noncoll)
                        {
                            states[st].psi[idx + pbasis] = 0.0;
                            states[st+1].psi[idx] = 0.0;
                            states[st+1].psi[idx + pbasis] = states[st].psi[idx];
                        }
                        idx++;

                    }               /* end for */
                }                   /* end for */
            }                       /* end for */

        }                           /* end for */

    }


}

