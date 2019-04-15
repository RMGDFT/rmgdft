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



#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "GpuAlloc.h"
#include "rmg_error.h"
#include "ErrorFuncs.h"

/*Call to this function needs to be preceeded by get_QI, since we use pct.Qidxptrlen,
 * which is setup in that function*/

static void reset_pct_arrays (int ion);

template void GetNlop<double> (Kpoint<double> **);
template void GetNlop<std::complex<double> > (Kpoint<std::complex<double>> **);

template <typename KpointType>
void GetNlop (Kpoint<KpointType> **Kptr)
{

    reset_pct_arrays (ct.num_ions);
    int projector_type = DELOCALIZED;
    if(ct.localize_projectors) projector_type = LOCALIZED;
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
    {
       if(Kptr[kpt]->BetaProjector) delete Kptr[kpt]->BetaProjector;
       Kptr[kpt]->BetaProjector = new Projector<KpointType>(Kptr[kpt], projector_type, pct.grid_npes, ct.num_ions);
    }
pct.num_tot_proj = Kptr[0]->BetaProjector->num_tot_proj;

    int num_nonloc_ions = Kptr[0]->BetaProjector->get_num_nonloc_ions();

    std::string newpath;

    if(ct.nvme_weights)
    {
        if(ct.nvme_weight_fd != -1) close(ct.nvme_weight_fd);
        if(ct.nvme_Bweight_fd != -1) close(ct.nvme_Bweight_fd);

        newpath = ct.nvme_weights_path + std::string("rmg_weight") + std::to_string(pct.spinpe) +
                  std::to_string(pct.kstart) + std::to_string(pct.gridpe);
        ct.nvme_weight_fd = FileOpenAndCreate(newpath, O_RDWR|O_CREAT|O_TRUNC, (mode_t)0600);
        
        newpath = ct.nvme_weights_path + std::string("rmg_Bweight") + std::to_string(pct.spinpe) +
                  std::to_string(pct.kstart) + std::to_string(pct.gridpe);
        ct.nvme_Bweight_fd = FileOpenAndCreate(newpath, O_RDWR|O_CREAT|O_TRUNC, (mode_t)0600);
    }

    int P0_BASIS = get_P0_BASIS();

    /* Grab some memory for temporary storage */
    int alloc = ct.max_nlpoints;

    if (alloc <get_NX_GRID() * get_NY_GRID() * get_NZ_GRID())
        alloc =get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

//    pct.num_tot_proj = num_nonloc_ions * ct.max_nl;

    pct.weight_size = pct.num_tot_proj * P0_BASIS + 128;

#if GPU_ENABLED
    cudaError_t custat;
    // Managed memory is faster when gpu memory is not constrained but 
    // pinned memory works better when it is constrained.
    if(ct.pin_nonlocal_weights)
    {
        custat = cudaMallocHost((void **)&pct.weight , pct.weight_size * sizeof(double));
        RmgCudaError(__FILE__, __LINE__, custat, "Error: cudaMallocHost failed.\n");
    }
    else
    {
        pct.weight = (double *)GpuMallocManaged(pct.weight_size * sizeof(double));
        int device = -1;
        cudaGetDevice(&device);
        cudaMemAdvise ( pct.weight, pct.weight_size * sizeof(double), cudaMemAdviseSetReadMostly, device);
    }
    for(size_t idx = 0;idx < pct.weight_size;idx++) pct.weight[idx] = 0.0;

    if(ct.need_Bweight) 
    {
        if(ct.pin_nonlocal_weights)
        {
            custat = cudaMallocHost((void **)&pct.Bweight , pct.weight_size * sizeof(double));
            RmgCudaError(__FILE__, __LINE__, custat, "Error: cudaMallocHost failed.\n");
        }
        else
        {
            pct.Bweight = (double *)GpuMallocManaged(pct.weight_size * sizeof(double));
            int device = -1;
            cudaGetDevice(&device);
            cudaMemAdvise ( pct.Bweight, pct.weight_size * sizeof(double), cudaMemAdviseSetReadMostly, device);
        }
        for(int idx = 0;idx < pct.weight_size;idx++) pct.Bweight[idx] = 0.0;
    }
    else {
        pct.Bweight = pct.weight;
    }
#else
    if(ct.nvme_weights)
    {
        pct.weight = (double *)CreateMmapArray(ct.nvme_weight_fd, pct.weight_size*sizeof(double));
        if(!pct.weight) rmg_error_handler(__FILE__,__LINE__,"Error: CreateMmapArray failed for weights. \n");
        madvise(pct.weight, pct.weight_size*sizeof(double), MADV_SEQUENTIAL);

        if(ct.need_Bweight) {
            pct.Bweight = (double *)CreateMmapArray(ct.nvme_Bweight_fd, pct.weight_size*sizeof(double));
            if(!pct.Bweight) rmg_error_handler(__FILE__,__LINE__,"Error: CreateMmapArray failed for bweights. \n");
        }
        else {
            pct.Bweight = pct.weight;
        }
    }
    else
    {
        pct.weight = new double[pct.weight_size]();
        if(ct.need_Bweight) {
            pct.Bweight = new double[pct.weight_size]();
        }
        else {
            pct.Bweight = pct.weight;
        }
    }
#endif


  
    
    // Set storage sequentially for real and imaginary components so we can transform storage pattern
#if GPU_ENABLED
    if (pct.newsintR_local)
        GpuFreeManaged(pct.newsintR_local);
#else
    if (pct.newsintR_local)
        delete [] pct.newsintR_local;
#endif
   
    int factor = 2;
    if(ct.is_gamma) factor = 1; 
    size_t sint_alloc = (size_t)(factor * ct.num_kpts_pe * num_nonloc_ions * ct.max_nl);
    sint_alloc *= (size_t)ct.max_states;
    sint_alloc += 16;    // In case of lots of vacuum make sure something is allocated otherwise allocation routine may fail
#if GPU_ENABLED
    pct.newsintR_local = (double *)GpuMallocManaged(sint_alloc * sizeof(double));
#else
    pct.newsintR_local = new double[sint_alloc]();
#endif

    KpointType *tsintnew_ptr = (KpointType *)pct.newsintR_local;

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++){
        Kptr[kpt]->sint_size = (size_t)num_nonloc_ions * (size_t)ct.max_states * (size_t)ct.max_nl;
        Kptr[kpt]->newsint_local = tsintnew_ptr;
        tsintnew_ptr += Kptr[kpt]->sint_size;
    }
    


} 


static void reset_pct_arrays (int num_ions)
{

    if (pct.weight != NULL) {
#if GPU_ENABLED
        if(ct.pin_nonlocal_weights)
        {
            cudaFreeHost(pct.weight);
        }
        else
        {
            cudaFree(pct.weight);
        }
#else
        if(ct.nvme_weights)
        {
            munmap(pct.weight, pct.weight_size*sizeof(double));
        }
        else
        {
            delete [] pct.weight;
        }
#endif
    }
    if ((pct.Bweight != NULL) && ct.need_Bweight) {
#if GPU_ENABLED
        if(ct.pin_nonlocal_weights)
        {
            cudaFreeHost(pct.Bweight);
        }
        else
        {
            cudaFree(pct.Bweight);
        }
#else
        if(ct.nvme_weights)
        {
            munmap(pct.Bweight, pct.weight_size*sizeof(double));
        }
        else
        {
            delete [] pct.Bweight;
        }
#endif
    }

}

