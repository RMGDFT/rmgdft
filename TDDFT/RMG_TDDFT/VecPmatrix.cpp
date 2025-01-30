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
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "blacs.h"
#include "RmgParallelFft.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

template void VecPHmatrix (Kpoint<double> *kptr, double *efield_tddft, int *desca, int tddft_start_state);
template void VecPHmatrix (Kpoint<std::complex<double>> *kptr, double *efield_tddft, int *desca, int tddft_start_state);

template <typename KpointType>
void VecPHmatrix (Kpoint<KpointType> *kptr, double *efield_tddft, int *desca, int tddft_start_state)
{

    int ictxt=desca[1], mb=desca[4], nb=desca[5], mxllda = desca[8];
    int mycol, myrow, nprow, npcol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    int num_states = kptr->nstates - tddft_start_state;
    int pbasis = kptr->pbasis;
    int pbasis_noncol = pbasis * ct.noncoll_factor;

    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));
    KpointType alpha(vel);
    KpointType beta(0.0);

    KpointType *block_matrix;

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
         trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }   

    //block_size = num_states;
    int num_blocks = (num_states + nb -1)/nb;
    // First time through allocate pinned memory for global_matrix1

    // 3 block matrix for px, py, pz operators
    int retval1 = MPI_Alloc_mem(3*num_states * nb * sizeof(KpointType) , MPI_INFO_NULL, &block_matrix);
    KpointType *block_matrix_x = block_matrix;
    KpointType *block_matrix_y = block_matrix_x + num_states * nb;
    KpointType *block_matrix_z = block_matrix_y + num_states * nb;


    if(retval1 != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
    }


    // V|psi> is in tmp_arrayT
    KpointType *psi = kptr->orbital_storage + tddft_start_state * pbasis_noncol;
    KpointType *psi_dev;
    if(kptr->psi_dev)
    {
        psi_dev = kptr->psi_dev + tddft_start_state * pbasis_noncol;
    }
    else
    {
        psi_dev = psi;
    }

    KpointType *psi_x = &kptr->orbital_storage[kptr->nstates * pbasis_noncol];  // use the memory of psi extra 3* state_block_size.
    KpointType *psi_y = &kptr->orbital_storage[kptr->nstates * pbasis_noncol] + nb * pbasis_noncol;  // use the memory of psi extra 3* state_block_size.
    KpointType *psi_z = &kptr->orbital_storage[kptr->nstates * pbasis_noncol] + 2*nb * pbasis_noncol;  // use the memory of psi extra 3* state_block_size.

    for(int ib = 0; ib < num_blocks; ib++)
    {
        int size_col = std::min(nb, num_states - ib * nb);
        int size_row = num_states - ib * nb;
        int this_block_size = size_col;


        for (int st1 = 0; st1 < size_col; st1++)
        {
            KpointType *psi1 = psi + (ib*nb + st1) * pbasis_noncol;
            KpointType *psi1_x = psi_x + st1 * pbasis_noncol;
            KpointType *psi1_y = psi_y + st1 * pbasis_noncol;
            KpointType *psi1_z = psi_z + st1 * pbasis_noncol;
            ApplyGradient(psi1, psi1_x, psi1_y, psi1_z, ct.force_grad_order, "Coarse");
            if(ct.noncoll)
            {
                ApplyGradient(psi1+pbasis, psi1_x+pbasis, psi1_y+pbasis, psi1_z+pbasis, ct.force_grad_order, "Coarse");
            }

            if(!ct.is_gamma)
            {
                std::complex<double> I_t(0.0, 1.0);
                std::complex<double> *psi_C, *psi_xC, *psi_yC, *psi_zC;
                psi_C = (std::complex<double> *) psi1;
                psi_xC = (std::complex<double> *) psi1_x;
                psi_yC = (std::complex<double> *) psi1_y;
                psi_zC = (std::complex<double> *) psi1_z;
                for(int i = 0; i < pbasis_noncol; i++)
                {
                    psi_xC[i] += I_t *  kptr->kp.kvec[0] * psi_C[i];
                    psi_yC[i] += I_t *  kptr->kp.kvec[1] * psi_C[i];
                    psi_zC[i] += I_t *  kptr->kp.kvec[2] * psi_C[i];
                }
            }

        }

        RmgGemm(trans_a, trans_n, size_row, size_col,  pbasis_noncol, alpha, psi_dev+ib*nb*pbasis_noncol, pbasis_noncol, psi_x, 
                pbasis_noncol, beta, block_matrix_x, size_row);
        BlockAllreduce((double *)block_matrix_x, (size_t)size_row * (size_t)size_col * (size_t)factor , pct.grid_comm);

        RmgGemm(trans_a, trans_n, size_row, size_col,  pbasis_noncol, alpha, psi_dev+ib*nb*pbasis_noncol, pbasis_noncol, psi_y, 
                pbasis_noncol, beta, block_matrix_y, size_row);
        BlockAllreduce((double *)block_matrix_y, (size_t)size_row * (size_t)size_col * (size_t)factor , pct.grid_comm);

        RmgGemm(trans_a, trans_n, size_row, size_col,  pbasis_noncol, alpha, psi_dev+ib*nb*pbasis_noncol, pbasis_noncol, psi_z, 
                pbasis_noncol, beta, block_matrix_z, size_row);
        BlockAllreduce((double *)block_matrix_z, (size_t)size_row * (size_t)size_col * (size_t)factor , pct.grid_comm);

        //block_matrix to distHij;
        if(myrow == ib%nprow)
        {
            int istart = (ib/nprow) *nb;
            for(int jb = ib; jb < num_blocks; jb++)
            {
                if(mycol == jb%npcol)
                    //  block (ib,jb) in this processor
                {
                    int this_block_size_col  = std::min(mb, num_states - mb * jb);
                    int jstart = (jb/npcol) * mb;
                    for(int i = 0; i < this_block_size; i++)
                    {
                        for(int j = 0; j < this_block_size_col; j++)
                        {
                            kptr->Pxmatrix_cpu[(jstart + j) * mxllda + i + istart] = block_matrix_x[ (j + (jb-ib) * mb ) * this_block_size + i];
                            kptr->Pymatrix_cpu[(jstart + j) * mxllda + i + istart] = block_matrix_y[ (j + (jb-ib) * mb ) * this_block_size + i];
                            kptr->Pzmatrix_cpu[(jstart + j) * mxllda + i + istart] = block_matrix_z[ (j + (jb-ib) * mb ) * this_block_size + i];
                        }
                    }
                }
            }
        }

        if(mycol == ib%npcol)
        {
            int istart = (ib/npcol) *nb;
            for(int jb = ib; jb < num_blocks; jb++)
            {
                if(myrow == jb%nprow)
                    //  block (jb,ib) in this processor
                {
                    int this_block_size_col  = std::min(mb, num_states - mb * jb);
                    int jstart = (jb/nprow) * mb;
                    for(int i = 0; i < this_block_size; i++)
                    {
                        for(int j = 0; j < this_block_size_col; j++)
                        {
                            kptr->Pxmatrix_cpu[(istart + i) * mxllda + j + jstart] = MyConj(block_matrix_x[ (j + (jb-ib) * mb ) * this_block_size + i]);
                            kptr->Pymatrix_cpu[(istart + i) * mxllda + j + jstart] = MyConj(block_matrix_y[ (j + (jb-ib) * mb ) * this_block_size + i]);
                            kptr->Pzmatrix_cpu[(istart + i) * mxllda + j + jstart] = MyConj(block_matrix_z[ (j + (jb-ib) * mb ) * this_block_size + i]);
                        }
                    }
                }

            }
        }

    }

    delete [] block_matrix;
}

