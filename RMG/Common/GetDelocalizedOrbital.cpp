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

    Projector<KpointType> *P = OrbitalProjector;
    size_t stride = P->get_pstride();

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

    int num_allorbitals_ldaUion = 0;
    for (size_t ion = 0; ion < this->ldaU->num_ldaU_ions; ++ion)
    {
        /* Generate ion pointer */
        int ion_idx = this->ldaU->ldaU_ion_index[ion];
        ION &Atom = Atoms[ion_idx];

        /* Get species type */
        SPECIES &AtomType = Species[Atom.species];
        num_allorbitals_ldaUion += AtomType.num_orbitals;
    }

    KpointType *npsi = (KpointType *)RmgMallocHost(num_allorbitals_ldaUion * pbasis * sizeof(KpointType));

    double coeff = 1.0;
    int wave_idx = 0;
    for (size_t ion = 0; ion < this->ldaU->num_ldaU_ions; ++ion)
    {
        /* Generate ion pointer */
        int ion_idx = this->ldaU->ldaU_ion_index[ion];
        ION &Atom = Atoms[ion_idx];

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
            std::cout << "LDA+U not working with localized orbitals" << std::endl;
            rmg_error_handler(__FILE__,__LINE__,"Used delocalized orbital for LDA+U, Terminating.");
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

    KpointType *weight;

    for (int ion = 0; ion < this->ldaU->num_ldaU_ions; ++ion)
    {
        int ion_idx = this->ldaU->ldaU_ion_index[ion];
        ION &Atom = Atoms[ion_idx];
        SPECIES &AtomType = Species[Atom.species];
        size_t offset = (size_t)ion * stride * (size_t)pbasis;
        weight = &orbital_weight[offset];
        //        std::complex<double> *Umm = AtomType.Umm;
        //        int lm_idx = 0;
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

    }

    RmgFreeHost(npsi);

    //  normalize the LDA+U orbitals

    double vel = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    vel = Rmg_L.get_omega() / vel;
    int num_orbitals = this->OrbitalProjector->get_num_tot_proj();
    if(ct.norm_conserving_pp)
    {
        for(int st = 0; st < num_orbitals; st++)
        {
            double sum = 0.0;
            for (int idx = 0; idx < pbasis; idx++)
            {
                sum += std::norm(orbital_weight[st * pbasis + idx] );
            }

            GlobalSums(&sum, 1, this->grid_comm);
            double tscale = vel * sum;
            // if the sum = 0, this orbital is not a LDA+U orbital and the values are zeros.
            if(tscale < 1.0e-5) continue;
            tscale = std::sqrt(1.0/tscale);

            if(ct.verbose && pct.imgpe == 0) std::cout << "norm constant " << tscale << std::endl;
            for (int idx = 0; idx < pbasis ; idx++)
            {
                orbital_weight[st * pbasis + idx] *= tscale;
            }

        }

    }
    else
    {

        size_t alloc = num_orbitals * pbasis;
#if HIP_ENABLED || CUDA_ENABLED
        this->BetaProjector->project(this, this->orbital_weight, this->newsint_local, 0,
                num_orbitals, this->nl_weight_gpu);
#else
        this->BetaProjector->project(this,this->orbital_weight, this->newsint_local, 0,
                num_orbitals, this->nl_weight);
#endif
        KpointType *tbuf = new KpointType[alloc]();
        AppS<KpointType>(this, this->newsint_local, this->orbital_weight, tbuf, 0, num_orbitals);

        for(int st = 0; st < num_orbitals; st++)
        {
            double sum = 0.0;
            for (int idx = 0; idx < pbasis; idx++)
            {
                sum += std::real(orbital_weight[st * pbasis + idx] * MyConj(tbuf[st * pbasis + idx]));
            }

            GlobalSums(&sum, 1, this->grid_comm);
            double tscale = vel * sum;
            if(tscale < 1.0e-5) continue;
            tscale = std::sqrt(1.0/tscale);

            if(ct.verbose && pct.imgpe == 0) std::cout << "norm constan for us " << tscale << std::endl;
            for (int idx = 0; idx < pbasis ; idx++)
            {
                orbital_weight[st * pbasis + idx] = tscale * tbuf[st * pbasis + idx];
            }

        }

        delete [] tbuf;
    }
}


