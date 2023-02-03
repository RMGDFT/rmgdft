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



template void HS_Scalapack (int nstates, int pbasis_noncoll, double *psi, double *hpsi, double *ns, int *desca, double *distHij, double *distSij);
template void HS_Scalapack (int nstates, int pbasis_noncoll, std::complex<double> *psi, std::complex<double> *hpsi, std::complex<double> *ns, int *desca, std::complex<double> *distHij, std::complex<double> *distSij);

template <typename KpointType>
void HS_Scalapack (int nstates, int pbasis_noncoll, KpointType *psi, KpointType *hpsi, KpointType *ns, int *desca, KpointType *distHij, KpointType *distSij)
{
    RmgTimer *RT1;

    int ictxt=desca[1], mb=desca[4], nb=desca[5], mxllda = desca[8];
    int mycol, myrow, nprow, npcol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);


    if( desca[2] != nstates || desca[3] != nstates)
    {
        printf("\n Warning: nstates %d != %d x %d matrix size \n", nstates, desca[2], desca[3]); 
    }
    if( mb != nb)
    {
        printf("\n Warning: block size difference %d != %d  \n", mb, nb);
    }
    // For Matrix multiplications
    double vel = Rmg_L.get_omega() /
        ((double)(Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1)));
    KpointType alphavel(vel);
    KpointType beta(0.0);

    // For MPI routines
    //int factor = 1;
    //if(!ct.is_gamma) factor = 2;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    if(typeid(KpointType) == typeid(std::complex<double>)) trans_a = trans_c;

    KpointType *block_matrix;
#if HIP_ENABLED || CUDA_ENABLED
    block_matrix = (KpointType *)GpuMallocHost( mb * nstates * sizeof(KpointType));
#else
    block_matrix = new KpointType[nstates * mb];
#endif


    //  For CPU only case and CUDA with managed memory psi_d is the same as orbital_storage but
    //  for HIP its a GPU buffer.
    KpointType *psi_dev = psi;
    KpointType *hpsi_dev = hpsi;
    KpointType *ns_dev = ns;
    if(ct.norm_conserving_pp)
    {
        ns_dev = psi_dev;
    }

#if HIP_ENABLED || CUDA_ENABLED
    if(ct.gpu_managed_memory == false && ct.use_cublasxt == false)
    {
        gpuMalloc((void **)&psi_dev, nstates * pbasis_noncoll * sizeof(KpointType));
        gpuMemcpy(psi_dev, psi, nstates * pbasis_noncoll * sizeof(KpointType), gpuMemcpyHostToDevice);
    }
#endif

    RT1 = new RmgTimer("4-Diagonalization: matrix");

    int num_blocks = (nstates + nb -1)/nb;

    for(int ib = 0; ib < num_blocks; ib++)
    {
        // upper triagle blocks only 
        // ib : (nstates - ib * nb) * block_size matrix
        // ib = 0: nstates * block_size matrix
        // ib = 1: (nstates - nb) * block_size matrix
        // block index (ib, ib:num_blocks)
        // this_block_size will be nb for the first num_blocks-1 block, the last block could be smaller than nb
        int this_block_size, length_block;

        length_block = nstates - nb * ib;
        this_block_size = std::min(nb, length_block);

        RmgTimer *RT1a = new RmgTimer("4-Diagonalization: matrix: Gemm");
        RmgGemm(trans_a, trans_n, this_block_size, length_block, pbasis_noncoll, alphavel, &hpsi_dev[ib*nb*pbasis_noncoll], pbasis_noncoll, 
                &psi_dev[ib*nb*pbasis_noncoll], pbasis_noncoll, beta, block_matrix, this_block_size);
        delete RT1a;

        RT1a = new RmgTimer("4-Diagonalization: matrix: Allreduce");
        BlockAllreduce(block_matrix, this_block_size * length_block , pct.grid_comm);
        delete RT1a;

        //block_matrix to distHij;
        if(myrow == ib%nprow)
        {
            int istart = (ib/nprow) *nb;
            for(int jb = ib; jb < num_blocks; jb++)
            {
                if(mycol == jb%npcol)
                    //  block (ib,jb) in this processor
                {
                    int this_block_size_col  = std::min(mb, nstates - mb * jb);
                    int jstart = (jb/npcol) * mb;
                    for(int i = 0; i < this_block_size; i++)
                    {
                        for(int j = 0; j < this_block_size_col; j++)
                        {
                            distHij[(jstart + j) * mxllda + i + istart] = block_matrix[ (j + (jb-ib) * mb ) * this_block_size + i];
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
                    int this_block_size_col  = std::min(mb, nstates - mb * jb);
                    int jstart = (jb/nprow) * mb;
                    for(int i = 0; i < this_block_size; i++)
                    {
                        for(int j = 0; j < this_block_size_col; j++)
                        {
                            distHij[(istart + i) * mxllda + j + jstart] = MyConj(block_matrix[ (j + (jb-ib) * mb ) * this_block_size + i]);
                        }
                    }
                }

            }
        }




        RT1a = new RmgTimer("4-Diagonalization: matrix: Gemm");

        RmgGemm(trans_a, trans_n, this_block_size, length_block, pbasis_noncoll, alphavel, &psi_dev[ib*nb*pbasis_noncoll], pbasis_noncoll, 
                &ns_dev[ib*nb*pbasis_noncoll], pbasis_noncoll, beta, block_matrix, this_block_size);
        delete RT1a;

        RT1a = new RmgTimer("4-Diagonalization: matrix: Allreduce");
        BlockAllreduce(block_matrix, this_block_size * length_block , pct.grid_comm);
        delete RT1a;

        //block_matrix to distHij;
        if(myrow == ib%nprow)
        {
            int istart = (ib/nprow) *nb;
            for(int jb = ib; jb < num_blocks; jb++)
            {
                if(mycol == jb%npcol)
                    //  block (ib,jb) in this processor
                {
                    int this_block_size_col  = std::min(mb, nstates - mb * jb);
                    int jstart = (jb/npcol) * mb;
                    for(int i = 0; i < this_block_size; i++)
                    {
                        for(int j = 0; j < this_block_size_col; j++)
                        {
                            distSij[(jstart + j) * mxllda + i + istart] = block_matrix[ (j + (jb-ib) * mb ) * this_block_size + i];
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
                    int this_block_size_col  = std::min(mb, nstates - mb * jb);
                    int jstart = (jb/nprow) * mb;
                    for(int i = 0; i < this_block_size; i++)
                    {
                        for(int j = 0; j < this_block_size_col; j++)
                        {
                            distSij[(istart + i) * mxllda + j + jstart] = MyConj(block_matrix[ (j + (jb-ib) * mb ) * this_block_size + i]);
                        }
                    }
                }

            }
        }

    }


    delete RT1;
#if HIP_ENABLED || CUDA_ENABLED
    if(ct.gpu_managed_memory == false && ct.use_cublasxt == false)
    {
        gpuFree(psi_dev);
        gpuFree(hpsi_dev);
    }
    GpuFreeHost(block_matrix);
#else
    delete [] block_matrix;
#endif

}
