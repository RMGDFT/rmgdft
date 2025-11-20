/*
 *
 * Copyright 2025 The RMG Project Developers. See the COPYRIGHT file 
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

#include "ortho.h"
#include "transition.h"
#include "RmgMatrix.h"

// Gram-Schmidt ortho for extra eigenvectors in Davidson solver.
// nbase: number of wavefunctions alreadt orthogonalized
// notcon: number of extra wavefunctions 
// psi: first nbase*this->pbasis will not be changed.
//      next notcon * this->pbasis will be orthogonalized. 
// only work for norm-conserving
// mat: matrix to hold <Psi|psi>, need to be allocated nbase * notconv at least
//


template ortho<double>::ortho(int, int);
template ortho<std::complex<double>>::ortho(int, int);
template ortho<double>::~ortho(void);
template ortho<std::complex<double>>::~ortho(void);
template void ortho<double>::orthogonalize(int, int, double *, bool);
template void ortho<std::complex<double>>::orthogonalize(int, int, std::complex<double>*, bool);


template <class T> ortho<T>::ortho(int max_states_in, int pbasis_in)
{
    if(!ct.norm_conserving_pp) return;
    this->max_states = max_states_in;
    this->pbasis = pbasis_in; 

#if HIP_ENABLED || CUDA_ENABLED
    size_t dfactor = 1;
    if(ct.kohn_sham_solver == DAVIDSON_SOLVER) dfactor = (size_t)ct.davidx;
    gpuMalloc((void **)&this->psi_d, dfactor * (size_t)this->max_states * 
              (size_t)this->pbasis * sizeof(T));
#endif
}

template <class T> ortho<T>::~ortho(void)
{
#if HIP_ENABLED || CUDA_ENABLED
    DeviceSynchronize();
    gpuFree(this->psi_d);
#endif
}

template <class T> void ortho<T>::orthogonalize(int nbase, int notcon, T *psi, bool dostage2)
{
    if(!ct.norm_conserving_pp) return;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    int factor = 1;
    if(typeid(T) == typeid(std::complex<double>)) {
        trans_a = trans_c;
        factor = 2;
    }

    size_t tlength = ((notcon + 2) * notcon / 2);
    size_t alloc = (notcon+nbase)*(notcon+nbase);
    T *mat = (T *)ct.get_gmatrix((alloc + tlength + 8192)*sizeof(T));
    size_t offset = 4096 * (alloc / 4096 + 1);
    T *tmat = mat + offset;

    double vel = Rmg_L.get_omega() /
        ((double)((size_t)Rmg_G->get_NX_GRID(1) * (size_t)Rmg_G->get_NY_GRID(1) * (size_t)Rmg_G->get_NZ_GRID(1)));
    T alphavel(vel);

    T zero(0.0);
    T one(1.0);
    T mone(-1.0);

#if HIP_ENABLED || CUDA_ENABLED
    T *psi_extra = &this->psi_d[nbase * this->pbasis];
    gpuMemcpy(this->psi_d, psi,
              (size_t)(notcon + nbase) * (size_t)this->pbasis * sizeof(T), gpuMemcpyHostToDevice);
#else
    T *psi_extra = &psi[nbase * this->pbasis];
    this->psi_d = psi;
#endif

    DeviceSynchronize();
    if(nbase > 0)
    {
        RmgTimer RT("MgridOrtho: 1st stage");
        // ortho to the first nbase states
        RmgGemm(trans_a, trans_n, nbase, notcon, this->pbasis, alphavel, this->psi_d, this->pbasis, psi_extra, this->pbasis, zero, mat, nbase);
        DeviceSynchronize();
        BlockAllreduce((double *)mat, (size_t)notcon*(size_t)nbase * (size_t)factor, pct.grid_comm);
        DeviceSynchronize();
        RmgGemm(trans_n, trans_n, this->pbasis, notcon, nbase, mone, this->psi_d, this->pbasis, mat, nbase, one, psi_extra, this->pbasis);
    }
    DeviceSynchronize();

    if(!dostage2)
    {
#if HIP_ENABLED || CUDA_ENABLED
        gpuMemcpy(&psi[nbase * this->pbasis], psi_extra,
              (size_t)notcon * (size_t)this->pbasis * sizeof(T), gpuMemcpyDeviceToHost);
#endif
        return;
    }

    RmgTimer *RT1 = new RmgTimer("MgridOrtho: 2nd stage overlaps");
#if HIP_ENABLED || CUDA_ENABLED
    T *mat_d;
    gpuMalloc((void **)&mat_d, (size_t)notcon * (size_t)notcon * sizeof(T));
#else   
    T *mat_d = mat;
#endif
    char *transt = "t";
    char *uplo = "u";
    char *diag = "N";
    char *side = "R";

    if constexpr (std::is_same_v<T, double>)
    {
        RmgSyrk( uplo, transt, notcon, this->pbasis, one, psi_extra, this->pbasis,
            zero, mat, notcon);
    }
    if constexpr (std::is_same_v<T, std::complex<double>>)
    {
        RmgGemm(trans_a, trans_n, notcon, notcon, this->pbasis, one, psi_extra,
                this->pbasis, psi_extra, this->pbasis, zero, mat, notcon);
    }
    delete RT1;

    /* get the global part */
    RT1 = new RmgTimer("MgridOrtho: allreduce");
    int length = factor * (notcon + 2) * notcon / 2;
    PackSqToTr("U", notcon, mat, notcon, tmat);
    BlockAllreduce((double *)tmat, length, pct.grid_comm);
    UnPackSqToTr("U", notcon, mat, notcon, tmat);
    delete RT1;

#if HIP_ENABLED || CUDA_ENABLED
    gpuMemcpy(mat_d, mat,
              (size_t)notcon * (size_t)notcon * sizeof(T), gpuMemcpyHostToDevice);
#endif

    /* compute the cholesky factor of the overlap matrix then subtract off projections */
    RT1 = new RmgTimer("MgridOrtho: cholesky");
    int info, info_trtri;
    T inv_vel(1.0/sqrt(vel));
    rmg_potrf(uplo, notcon, mat_d, notcon, &info);
    delete RT1;
    RT1 = new RmgTimer("MgridOrtho: inverse");
    rmg_trtri(uplo, diag, notcon, mat_d, notcon, &info_trtri);
    delete RT1;
    RT1 = new RmgTimer("MgridOrtho: 2nd stage update");
    rmg_trmm(side, uplo, "N", diag, this->pbasis, notcon, inv_vel, mat_d, notcon, psi_extra, this->pbasis);
    delete RT1;

#if HIP_ENABLED || CUDA_ENABLED
    gpuMemcpy(&psi[nbase * this->pbasis], psi_extra,
          (size_t)notcon * (size_t)this->pbasis * sizeof(T), gpuMemcpyDeviceToHost);
    gpuFree(mat_d);
#endif

    if (info != 0)
        throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << ". Matrix not positive definite or argument error. Terminating.\n";

    if (info_trtri != 0)
    throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << "info = " <<info << ". tritri problem.\n";

}


