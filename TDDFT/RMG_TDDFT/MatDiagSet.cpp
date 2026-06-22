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

// set the matrix with diagonal elements to be diag_elem, and others to 0.0
template void MatDiagSet<double> (double *mat, std::vector<double> diag_elem, double beta, int numst, Scalapack &SP);
template void MatDiagSet<std::complex<double>> (std::complex<double> *mat, std::vector<double> diag_elem, double beta, int numst, Scalapack &SP);
template <typename MatrixType>
    void MatDiagSet (MatrixType *mat,  std::vector<double> diag_elem, double beta, int numst, Scalapack &SP)
{

    if(ct.tddft_tiledMM)
    {
        int numst_pe = numst/pct.local_comm_npes;
        for(int i = 0; i < numst_pe; i++)
        {
            mat[i * numst + i+pct.local_rank * numst_pe] += beta * diag_elem[i+pct.local_rank * numst_pe];
        }
    }
    else
    {
        int *desca = SP.GetDistDesca();
        int ictxt=desca[1], mb=desca[4], nb=desca[5], mxllda = desca[8];
        int mycol, myrow, nprow, npcol;
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
        int izero = 0;
        for(int i = 0; i < SP.GetDistMdim(); i++)
        {
            int i1 = i+1;
            int i_glob = indxl2g(&i1, &mb, &myrow, &izero, &nprow);
            for(int j = 0; j < SP.GetDistNdim(); j++)
            {
                int j1 = j+1;
                int j_glob = indxl2g(&j1, &nb, &mycol, &izero, &npcol);
                if(i_glob == j_glob)
                {
                    mat[i + j * mxllda] += beta * diag_elem[i_glob-1];
                }
            }
        }

    }
}
// set the matrix with diagonal elements to be diag_elem, and others to 0.0
template void MatDiagGet<double> (double *mat, std::vector<double> diag_elem, int numst, Scalapack &SP);
template void MatDiagGet<std::complex<double>> (std::complex<double> *mat, std::vector<double> diag_elem, int numst, Scalapack &SP);
template <typename MatrixType>
    void MatDiagGet (MatrixType *mat,  std::vector<double> diag_elem, int numst, Scalapack &SP)
{

    std::fill(diag_elem.begin(), diag_elem.end(), 0.0);
    if(ct.tddft_tiledMM)
    {
        int numst_pe = numst/pct.local_comm_npes;
        for(int i = 0; i < numst_pe; i++)
        {
            diag_elem[i+pct.local_rank * numst_pe] = 
                std::real(mat[i * numst + i+pct.local_rank * numst_pe]);
        }

        rmg::all_reduce(diag_elem.data(), numst, pct.local_comm);

    }
    else
    {
        int *desca = SP.GetDistDesca();
        int ictxt=desca[1], mb=desca[4], nb=desca[5], mxllda = desca[8];
        int mycol, myrow, nprow, npcol;
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
        int izero = 0;
        for(int i = 0; i < SP.GetDistMdim(); i++)
        {
            int i1 = i+1;
            int i_glob = indxl2g(&i1, &mb, &myrow, &izero, &nprow);
            for(int j = 0; j < SP.GetDistNdim(); j++)
            {
                int j1 = j+1;
                int j_glob = indxl2g(&j1, &nb, &mycol, &izero, &npcol);
                if(i_glob == j_glob)
                {
                    diag_elem[i_glob-1] = std::real(mat[i + j * mxllda]);
                }
            }
        }

        rmg::all_reduce(diag_elem.data(), numst, SP.GetComm());

    }
}
