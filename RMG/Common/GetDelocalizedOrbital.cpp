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
#include "GlobalSums.h"
#include "params.h"
#include "rmg_error.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"


// Used to generate LDA+U orbital projectors that span the full space
template void Kpoint<double>::GetDelocalizedOrbital(void);
template void Kpoint<std::complex<double>>::GetDelocalizedOrbital(void);
template <class KpointType> void Kpoint<KpointType>::GetDelocalizedOrbital (void)
{

    
    double *kvec = kp.kvec;


    for(size_t idx = 0; idx < orbital_weight_size; idx++)
            orbital_weight[idx] = 0.0;


    /* Loop over ions */

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

    int state_count = CountAtomicOrbitals();

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

    wave_idx = 0;

    Projector<KpointType> *P = OrbitalProjector;
    size_t stride = P->get_pstride();
    KpointType *weight;

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        ION &Atom = Atoms[ion];
        SPECIES &AtomType = Species[Atom.species];
        if(!AtomType.is_ldaU) 
        {
            wave_idx += AtomType.num_orbitals;
            continue;
        }
        size_t offset = (size_t)ion * stride * (size_t)pbasis;
        weight = &orbital_weight[offset];
        std::complex<double> *Umm = AtomType.Umm;
        int lm_idx = 0;
        for (int ip = 0; ip < AtomType.num_atomic_waves; ip++)
        {
            int li = AtomType.atomic_wave_l[ip];
            double ji = AtomType.atomic_wave_j[ip];

            if(strcasecmp(AtomType.ldaU_label.c_str(), AtomType.atomic_wave_label[ip].c_str() ))
            {
                wave_idx += 2*li+1;
                continue;
            }
            if(ct.spinorbit || ct.noncoll)
            {
                double factor_j = 1.0;
                if(AtomType.is_spinorb)
                    factor_j = (2.0 * ji + 1.0) /(2.0 * li +1 );

                for(int m = 0; m < 2*li+1; m++)
                {
                    std::complex<double> *psi = (std::complex<double> *)&weight[m * pbasis];

                    for(int idx = 0; idx < pbasis; idx++)
                    {
                        psi[idx] += factor_j * npsi[(wave_idx+m) * pbasis + idx];
                    }
                }

                wave_idx += 2 * li + 1;
            }

            else 
            {
                for (int m=0; m < 2*li+1; m++)
                {
                    for (int idx = 0; idx < pbasis; idx++) 
                        weight[m * pbasis + idx] = npsi[wave_idx * pbasis + idx];;
                    wave_idx ++;

                }
            }
        }
        
// if spin orbit, normalize the averaged orbitals
        double vel = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
        vel = Rmg_L.get_omega() / vel;
        if(AtomType.is_spinorb)
        {
            for(int st = 0; st < stride; st++)
            {
                double sum = 0.0;
                for (int idx = 0; idx < pbasis; idx++)
                {
                    sum += std::norm(weight[st * pbasis + idx] );
                }

                GlobalSums(&sum, 1, this->grid_comm);
                double tscale = vel * sum;
                tscale = std::sqrt(1.0/tscale);

                for (int idx = 0; idx < pbasis ; idx++)
                {
                    weight[st * pbasis + idx] *= tscale;
                }

            }
        }
    }
    GpuFreeManaged(npsi);
}


