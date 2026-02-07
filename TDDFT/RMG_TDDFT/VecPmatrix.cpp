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
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "rmg_reduce.h"
#include "Kpoint.h"
#include "rmg_gemm.h"
#include "Subdiag.h"
#include "GpuAlloc.h"

#include "blas.h"
#include "blacs.h"
#include "RmgParallelFft.h"
#include "RmgException.h"


#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

template void VecPHmatrix<double> (Kpoint<double> *kptr, double *efield_tddft, int *desca, int tddft_start_state, int num_states);
template void VecPHmatrix<std::complex<double>> (Kpoint<std::complex<double>> *kptr, double *efield_tddft, int *desca, int tddft_start_state, int num_states);
template <typename OrbitalType>
    void VecPHmatrix (Kpoint<OrbitalType> *kptr, double *efield_tddft, int *desca, int tddft_start_state, int num_states)
{

    int ictxt=desca[1], mb=desca[4], nb=desca[5], mxllda = desca[8];
    int mycol, myrow, nprow, npcol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    rmg::grid *G = kptr->G;
    Lattice *L = kptr->L;

    int pbasis = kptr->pbasis;
    int pbasis_noncol = pbasis * ct.noncoll_factor;

    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));
    //  alpha take care of i in moment operator
    OrbitalType alpha(vel);
    OrbitalType beta(0.0);
    std::complex<double> I_t(0.0, 1.0);

    OrbitalType *block_matrix;

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    trans_a = trans_c;

    //block_size = num_states;
    int num_blocks = (num_states + nb -1)/nb;
    // First time through allocate pinned memory for global_matrix1

    // 3 block matrix for px, py, pz operators
    int retval1 = MPI_Alloc_mem(3*num_states * nb * sizeof(OrbitalType) , MPI_INFO_NULL, &block_matrix);
    OrbitalType *block_matrix_x = block_matrix;
    OrbitalType *block_matrix_y = block_matrix_x + num_states * nb;
    OrbitalType *block_matrix_z = block_matrix_y + num_states * nb;

    if(retval1 != MPI_SUCCESS) {
        rmg::error("Memory allocation failure in HmatrixUpdate");
    }


    // V|psi> is in tmp_arrayT
    OrbitalType *psi = kptr->orbital_storage + tddft_start_state * pbasis_noncol;
    OrbitalType *psi_dev;
    if(kptr->psi_dev)
    {
        psi_dev = kptr->psi_dev + tddft_start_state * pbasis_noncol;
    }
    else
    {
        psi_dev = psi;
    }

    OrbitalType *psi_x = &kptr->orbital_storage[kptr->nstates * pbasis_noncol];  // use the memory of psi extra 3* state_block_size.
    OrbitalType *psi_y = &kptr->orbital_storage[kptr->nstates * pbasis_noncol] + nb * pbasis_noncol;  // use the memory of psi extra 3* state_block_size.
    OrbitalType *psi_z = &kptr->orbital_storage[kptr->nstates * pbasis_noncol] + 2*nb * pbasis_noncol;  // use the memory of psi extra 3* state_block_size.

    for(int ib = 0; ib < num_blocks; ib++)
    {
        // upper triagle blocks only
        // ib : (nstates - ib * nb) * block_size matrix
        // ib = 0: nstates * block_size matrix
        // ib = 1: (nstates - nb) * block_size matrix
        // block index (ib, ib:num_blocks)
        // this_block_size will be nb for the first num_blocks-1 block, the last block could be smaller than nb


        int this_block_size, length_block;
        length_block = num_states - nb * ib;
        this_block_size = std::min(nb, length_block);


        for (int st1 = 0; st1 < this_block_size; st1++)
        {
            OrbitalType *psi1 = psi + (ib*nb + st1) * pbasis_noncol;
            OrbitalType *psi1_x = psi_x + st1 * pbasis_noncol;
            OrbitalType *psi1_y = psi_y + st1 * pbasis_noncol;
            OrbitalType *psi1_z = psi_z + st1 * pbasis_noncol;
            ApplyGradient(psi1, psi1_x, psi1_y, psi1_z, ct.force_grad_order, "Coarse");
            if(ct.noncoll)
            {
                ApplyGradient(psi1+pbasis, psi1_x+pbasis, psi1_y+pbasis, psi1_z+pbasis, ct.force_grad_order, "Coarse");
            }

            if(!ct.is_gamma)
            {
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

        rmg::gemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_x, pbasis_noncol, psi_dev,
                pbasis_noncol, beta, block_matrix_x, this_block_size);
        rmg::block_reduce((double *)block_matrix_x, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

        rmg::gemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_y, pbasis_noncol, psi_dev,
                pbasis_noncol, beta, block_matrix_y, this_block_size);
        rmg::block_reduce((double *)block_matrix_y, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

        rmg::gemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_z, pbasis_noncol, psi_dev,
                pbasis_noncol, beta, block_matrix_z, this_block_size);
        rmg::block_reduce((double *)block_matrix_z, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

        //block_matrix to distHij;
        if(ct.tddft_tiledMM)
        {
            int istart = ib *nb;
            for(int i = 0; i < this_block_size; i++)
            {
                for(int j = 0; j < num_states/pct.local_comm_npes; j++)
                {
                    int jglob = j + pct.local_rank * num_states/pct.local_comm_npes ;
                    kptr->Pxmatrix_cpu[ j * num_states + i + istart] = I_t*block_matrix_x[ jglob * this_block_size + i];
                    kptr->Pymatrix_cpu[ j * num_states + i + istart] = I_t*block_matrix_y[ jglob * this_block_size + i];
                    kptr->Pzmatrix_cpu[ j * num_states + i + istart] = I_t*block_matrix_z[ jglob * this_block_size + i];
                }
            }

        }
        else
        {
            if(myrow == ib%nprow)
            {
                int istart = (ib/nprow) *nb;
                for(int jb = 0; jb < num_blocks; jb++)
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
                                kptr->Pxmatrix_cpu[(jstart + j) * mxllda + i + istart] = I_t*block_matrix_x[ (j + jb * mb ) * this_block_size + i];
                                kptr->Pymatrix_cpu[(jstart + j) * mxllda + i + istart] = I_t*block_matrix_y[ (j + jb * mb ) * this_block_size + i];
                                kptr->Pzmatrix_cpu[(jstart + j) * mxllda + i + istart] = I_t*block_matrix_z[ (j + jb * mb ) * this_block_size + i];
                            }
                        }
                    }
                }
            }
        }

    }

    delete [] block_matrix;
    //  if(pct.worldrank == 4)
    //  {
    //      printf("\n %d %d ", myrow, mycol);

    //      printf("\nkvec %f %f %f\n", kptr->kp.kvec[0] , kptr->kp.kvec[1] , kptr->kp.kvec[2] );
    //      fflush(NULL);
    //      for(int i = 0; i < 8; i++)
    //      {
    //          printf("\n aaa ");
    //          for(int j = 0; j < 8; j++)
    //              printf(" %8.3e ", std::real(kptr->Pxmatrix_cpu[i *8 + j]));
    //      }
    //      fflush(NULL);
    //      for(int i = 0; i < 8; i++)
    //      {
    //          printf("\n bbb ");
    //          for(int j = 0; j < 8; j++)
    //              printf(" %8.3e ", std::imag(kptr->Pxmatrix_cpu[i *8 + j]));
    //      }
    //      fflush(NULL);
    //  }

}

