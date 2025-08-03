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
#include "RmgException.h"


#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

void CurrentNlpp (Kpoint<double> *kptr, int *desca, int tddft_start_state)
{
    throw RmgFatalException() << "TDDFT vector potential mode: wave function cannot be real now " << "\n";
}

void CurrentNlpp (Kpoint<std::complex<double>> *kptr, int *desca, int tddft_start_state)
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
    //  alpha take care of i in moment operator
    std::complex<double> alpha(0.0, vel);
    std::complex<double> beta(0.0);

    std::complex<double> *block_matrix;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    trans_a = trans_c;

    //block_size = num_states;
    int num_blocks = (num_states + nb -1)/nb;
    // First time through allocate pinned memory for global_matrix1

    // 3 block matrix for px, py, pz operators
    int retval1 = MPI_Alloc_mem(3*num_states * nb * sizeof(std::complex<double>) , MPI_INFO_NULL, &block_matrix);
    std::complex<double> *block_matrix_x = block_matrix;
    std::complex<double> *block_matrix_y = block_matrix_x + num_states * nb;
    std::complex<double> *block_matrix_z = block_matrix_y + num_states * nb;

    if(nb != ct.scalapack_block_factor)
    {
        rmg_error_handler (__FILE__, __LINE__, "state_block_size must be same as scalack_block_factor\n");
    }
    if(retval1 != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
    }


    // V|psi> is in tmp_arrayT
    std::complex<double> *psi = kptr->orbital_storage + tddft_start_state * pbasis_noncol;
    std::complex<double> *psi_dev;
    if(kptr->psi_dev)
    {
        psi_dev = kptr->psi_dev + tddft_start_state * pbasis_noncol;
    }
    else
    {
        psi_dev = psi;
    }

    std::complex<double> *psi_x = &kptr->orbital_storage[kptr->nstates * pbasis_noncol];  // use the memory of psi extra 3* state_block_size.
    std::complex<double> *psi_y = &kptr->orbital_storage[kptr->nstates * pbasis_noncol] + nb * pbasis_noncol;  // use the memory of psi extra 3* state_block_size.
    std::complex<double> *psi_z = &kptr->orbital_storage[kptr->nstates * pbasis_noncol] + 2*nb * pbasis_noncol;  // use the memory of psi extra 3* state_block_size.

    std::complex<double> *ns = kptr->ns;
    std::complex<double> *nv = kptr->nv;
    std::complex<double> *newsint_local = kptr->newsint_local;


    int factor = 2;
    int ix, iy, iz;
    Rmg_G->pe2xyz (pct.gridpe, &ix, &iy, &iz);
    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);

    int px0_grid = Rmg_G->get_PX0_GRID(1);
    int py0_grid = Rmg_G->get_PY0_GRID(1);
    int pz0_grid = Rmg_G->get_PZ0_GRID(1);
    double xoff = ix * px0_grid * hxgrid;
    double yoff = iy * py0_grid * hygrid;
    double zoff = iz * pz0_grid * hzgrid;

    double xtal[3], xcrt[3];

    for(int ib = 0; ib < num_blocks; ib++)
    {
        // upper triagle blocks only
        // ib : (nstates - ib * nb) * block_size matrix
        // ib = 0: nstates * block_size matrix
        // ib = 1: (nstates - nb) * block_size matrix
        // block index (ib, ib:num_blocks)
        // this_block_size will be nb for the first num_blocks-1 block, the last block could be smaller than nb

        int this_block_size, length_block;
        length_block = num_states - ib * nb;
        this_block_size = std::min(nb, length_block);
        int st_start = ib *nb + tddft_start_state;
        int st_end   = st_start + this_block_size;

        AppNls(kptr, newsint_local, kptr->Kstates[0].psi, nv, ns, st_start, st_end);

        for (int st1 = 0; st1 < this_block_size; st1++)
        {

            for(int i = 0; i < px0_grid; i++)
            {
                for(int j = 0; j < py0_grid; j++)
                {
                    for(int k = 0; k < pz0_grid; k++)
                    {

                        xtal[0] = xoff + i * hxgrid;
                        xtal[1] = yoff + j * hygrid;
                        xtal[2] = zoff + k * hzgrid;
                        Rmg_L.to_cartesian(xtal, xcrt);

                        int idx = st1 * pbasis_noncol + i * py0_grid * pz0_grid + j * pz0_grid + k;
                        psi_x[idx] = nv[idx] * xcrt[0];
                        psi_y[idx] = nv[idx] * xcrt[1];
                        psi_z[idx] = nv[idx] * xcrt[2];

                    } 
                }
            }


            if(ct.noncoll)
            {
                rmg_printf("\nAAAAAA\n"); 
                exit(0);
            }

        }

        RmgGemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_x, pbasis_noncol, psi_dev,
                pbasis_noncol, beta, block_matrix_x, this_block_size);
        BlockAllreduce((double *)block_matrix_x, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

        RmgGemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_y, pbasis_noncol, psi_dev,
                pbasis_noncol, beta, block_matrix_y, this_block_size);
        BlockAllreduce((double *)block_matrix_y, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

        RmgGemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_z, pbasis_noncol, psi_dev,
                pbasis_noncol, beta, block_matrix_z, this_block_size);
        BlockAllreduce((double *)block_matrix_z, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

        //block_matrix to distHij;
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
                            kptr->Pxmatrix_cpu[(jstart + j) * mxllda + i + istart] += block_matrix_x[ (j + jb * mb ) * this_block_size + i];
                            kptr->Pymatrix_cpu[(jstart + j) * mxllda + i + istart] += block_matrix_y[ (j + jb * mb ) * this_block_size + i];
                            kptr->Pzmatrix_cpu[(jstart + j) * mxllda + i + istart] += block_matrix_z[ (j + jb * mb ) * this_block_size + i];
                        }
                    }
                }
            }
        }

        if(mycol == ib%npcol)
        {
            int istart = (ib/npcol) *nb;
            for(int jb = 0; jb < num_blocks; jb++)
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
                            kptr->Pxmatrix_cpu[(istart + i) * mxllda + j + jstart] += MyConj(block_matrix_x[ (j + jb * mb ) * this_block_size + i]);
                            kptr->Pymatrix_cpu[(istart + i) * mxllda + j + jstart] += MyConj(block_matrix_y[ (j + jb * mb ) * this_block_size + i]);
                            kptr->Pzmatrix_cpu[(istart + i) * mxllda + j + jstart] += MyConj(block_matrix_z[ (j + jb * mb ) * this_block_size + i]);
                        }
                    }
                }
            }
        }

    }

    delete [] block_matrix;
//  rmg_printf("kvec %f", kptr->kp.kvec[0] );
//  for(int i = 0; i < 8; i++)
//  {
//      rmg_printf("\n aaa ");
//      for(int j = 0; j < 8; j++)
//          rmg_printf(" %8.3e ", std::real(kptr->Pxmatrix_cpu[i *8 + j]));
//  }
//  for(int i = 0; i < 8; i++)
//  {
//      rmg_printf("\n bbb ");
//      for(int j = 0; j < 8; j++)
//          rmg_printf(" %8.3e ", std::imag(kptr->Pxmatrix_cpu[i *8 + j]));
//  }

}
