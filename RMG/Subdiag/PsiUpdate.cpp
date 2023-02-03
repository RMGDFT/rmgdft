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
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "GpuAlloc.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "blas.h"
#include "blacs.h"
#include "RmgMatrix.h"
#include "Functional.h"


#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"


template <typename KpointType>
void PsiUpdate (int nstates, int pbasis_noncoll, KpointType *distAij, int *desca, KpointType *psi, KpointType *hpsi, KpointType *matrix_diag);

template void PsiUpdate<double> (int nstates, int pbasis_noncoll, double *distAij, int *desca, 
        double *psi, double *hpsi, double *matrix_diag);
template void PsiUpdate<std::complex<double>> (int nstates, int pbasis_noncoll, std::complex<double> *distAij, int *desca, 
        std::complex<double> *psi, std::complex<double> *hpsi, std::complex<double> *matrix_diag);

template <typename KpointType>
void PsiUpdate (int nstates, int pbasis_noncoll, KpointType *distAij, int *desca, KpointType *psi, KpointType *hpsi, KpointType *matrix_diag)
{
    RmgTimer *RT1;

    int ictxt=desca[1], mb=desca[4], nb=desca[5], mxllda = desca[8];
    int mycol, myrow, nprow, npcol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    KpointType alpha(1.0);
    KpointType beta(0.0);
    // Begin rotation
    //
    char *trans_n = "n";
    KpointType *block_matrix;
    KpointType *block_matrix_dev;
    KpointType *psi_dev;
    KpointType *psi_new_dev;


#if HIP_ENABLED || CUDA_ENABLED
    block_matrix = (KpointType *)GpuMallocHost( nb * nstates * sizeof(KpointType));
    gpuMalloc((void **)&psi_new_dev, nb * pbasis_noncoll * sizeof(KpointType));
#else
    block_matrix = new KpointType[nstates * nb];
    psi_new_dev = new KpointType[nb*pbasis_noncoll];
#endif

    block_matrix_dev = block_matrix;
    psi_dev = psi;
    

    int num_blocks = (nstates + nb -1)/nb;

    for(int ib = 0; ib < num_blocks; ib++)
    {
        int this_block_size_row;
        this_block_size_row = std::min(mb, nstates - mb * ib);

        RT1 = new RmgTimer("4-Diagonalization: Update orbitals: gather");
        size_t size_mat = this_block_size_row * nstates;
        for(size_t i = 0; i < size_mat; i++) block_matrix[i] = 0.0;
        for(int jb = 0; jb < num_blocks; jb++)
        {
            int this_block_size_col;
            this_block_size_col = std::min(nb, nstates - nb * jb);

            if(myrow == jb%nprow && mycol == ib%npcol)
            {
                for(int i = 0; i < this_block_size_row; i++)
                {
                    for(int j = 0; j < this_block_size_col; j++)
                    {
                        block_matrix[i * nstates + jb * nb + j] = distAij[ ((ib/npcol) * mb + i) * mxllda + (jb/nprow) * nb + j];

                    }
                }
            }
        }

        BlockAllreduce(block_matrix, size_mat, pct.grid_comm);
//        if(pct.imgpe == 0 && pct.gridpe == 0) 
//            for(int i = 0; i < this_block_size_row; i++)
//                for(int j = 0; j < nstates; j++) printf("\n %d %d %f eee", i+ib*nb, j,block_matrix[i*nstates + j]); 

        for(int i = 0; i < this_block_size_row; i++) matrix_diag[i] = block_matrix[i * nstates + ib*mb + i];
        delete RT1;



        RT1 = new RmgTimer("4-Diagonalization: Update orbitals: gemm");
        RmgGemm(trans_n, trans_n, pbasis_noncoll, this_block_size_row, nstates, alpha, 
                psi_dev, pbasis_noncoll, block_matrix_dev, nstates, 
                beta, psi_new_dev, pbasis_noncoll);
        delete RT1;

#if HIP_ENABLED || CUDA_ENABLED
    gpuMemcpy(&hpsi[ib*nb*pbasis_noncoll], psi_new_dev, this_block_size_row * pbasis_noncoll * sizeof(KpointType), gpuMemcpyDeviceToHost);
#else
    memcpy(&hpsi[ib*nb*pbasis_noncoll], psi_new_dev, this_block_size_row * pbasis_noncoll * sizeof(KpointType));

#endif
                
    }

#if HIP_ENABLED || CUDA_ENABLED
    GpuFreeHost(block_matrix);
    gpuFree(psi_new_dev);
#else
    delete [] block_matrix;
    delete [] psi_new_dev;
#endif

}
