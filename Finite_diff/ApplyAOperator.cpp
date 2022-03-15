/*
 *
 * Copyright (c) 2014, Emil Briggs
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
#include "const.h"
#include "TradeImages.h"
#include "RmgException.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "transition.h"
#include "RmgParallelFft.h"
#include "GpuAlloc.h"
#include "Gpufuncs.h"
#include "ErrorFuncs.h"
#include "rmg_complex.h"

//  Applies the A operator to a wavefunction. The A operator is defined as
//
//  A(psi) = Laplacian(psi) + 2.0*i*[k dot Gradient(psi)]
//
//  For gamma this is just the laplacian.


template double ApplyAOperator<float>(float *, float *, int, int, int, double, double, double, int, double *kvec);
template double ApplyAOperator<double>(double *, double *,int, int, int, double, double, double, int, double *kvec);
template double ApplyAOperator<std::complex<float> >(std::complex<float> *, std::complex<float> *,int, int, int, double, double, double, int, double *kvec);
template double ApplyAOperator<std::complex<double> >(std::complex<double> *, std::complex<double> *,int, int, int, double, double, double, int, double *kvec);

template double ApplyAOperator<float>(float *, float *);
template double ApplyAOperator<double>(double *, double *);
template double ApplyAOperator<std::complex<float> >(std::complex<float> *, std::complex<float> *);
template double ApplyAOperator<std::complex<double> >(std::complex<double> *, std::complex<double> *);

template double ApplyAOperator<float>(float *, float *, double *);
template double ApplyAOperator<double>(double *, double *, double *);
template double ApplyAOperator<std::complex<float> >(std::complex<float> *, std::complex<float> *, double *);
template double ApplyAOperator<std::complex<double> >(std::complex<double> *, std::complex<double> *, double *);

template double ApplyAOperator<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int);
template double ApplyAOperator<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);


static void *rbufs[MAX_RMG_THREADS];

template <typename DataType>
double ApplyAOperator (DataType *a, DataType *b, double *kvec)
{
    int density = 1;

    int dimx = Rmg_G->get_PX0_GRID(density);
    int dimy = Rmg_G->get_PY0_GRID(density);
    int dimz = Rmg_G->get_PZ0_GRID(density);

    double gridhx = Rmg_G->get_hxgrid(density);
    double gridhy = Rmg_G->get_hygrid(density);
    double gridhz = Rmg_G->get_hzgrid(density);
                                                              
    return ApplyAOperator (a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, kvec);

}


template <typename DataType>
double ApplyAOperator (DataType *a, DataType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order, double *kvec)
{
    BaseThread *Th = BaseThread::getBaseThread(0);
    int tid = Th->get_thread_tid();
    if(tid < 0) tid = 0;


    if(ct.kohn_sham_ke_fft || Rmg_L.get_ibrav_type() == None)
    {
        FftLaplacianCoarse(a, b);    
        if(!ct.is_gamma)
        {
            int pbasis = dimx*dimy*dimz;
            DataType *gx = new DataType[pbasis];
            DataType *gy = new DataType[pbasis];
            DataType *gz = new DataType[pbasis];
            std::complex<double> *kdr = new std::complex<double>[pbasis];

            FftGradientCoarse(a, gx, gy, gz);
            // if complex orbitals compute dot products as well
            std::complex<double> I_t(0.0, 1.0);
            for(int idx = 0;idx < pbasis;idx++)
            {
                kdr[idx] = -I_t * (kvec[0] * std::complex<double>(std::real(gx[idx]), std::imag(gx[idx])) +
                                   kvec[1] * std::complex<double>(std::real(gy[idx]), std::imag(gy[idx])) +
                                   kvec[2] * std::complex<double>(std::real(gz[idx]), std::imag(gz[idx])));
            }
            if(typeid(DataType) == typeid(std::complex<double>))
            {
                std::complex<double> *pptr = (std::complex<double> *)b;
                for(int idx = 0;idx < pbasis;idx++) pptr[idx] = pptr[idx] - 2.0*kdr[idx];
            }
            if(typeid(DataType) == typeid(std::complex<float>))
            {
                std::complex<float> *pptr = (std::complex<float> *)b;
                for(int idx = 0;idx < pbasis;idx++) pptr[idx] = pptr[idx] - 2.0*kdr[idx];
            }
            delete [] kdr;
            delete [] gz;
            delete [] gy;
            delete [] gx;
        }

        FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
        DataType *ptr = NULL;
        // When ptr=NULL this does not do the finite differencing but just
        // returns the value of the diagonal element.
        double fd_diag = FD.app8_del2 (ptr, ptr, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        return 2.0*fd_diag;
    }

    double cc = 0.0;
    FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
    int sbasis = (pct.coalesce_factor*Rmg_G->get_PX0_GRID(1) + order) * (dimy + order) * (dimz + order);
    size_t alloc = (sbasis + 64) * sizeof(double);
    if(!ct.is_gamma) alloc *= 2;
    int images = order / 2;

    // Allocate per thread buffers if not done yet
    if(!rbufs[tid])
    {
#if HIP_ENABLED || CUDA_ENABLED
#if HIP_ENABLED
        hipSetDevice(ct.hip_dev);
#endif
#if CUDA_ENABLED
        cudaSetDevice(ct.cu_dev);
#endif
        gpuMallocHost((void **)&rbufs[tid], alloc);
#else
        MPI_Alloc_mem(alloc, MPI_INFO_NULL, &rbufs[tid]);
#endif
    }
    DataType *rptr = (DataType *)rbufs[tid];

    int special = ((Rmg_L.get_ibrav_type() == HEXAGONAL) ||
                   (Rmg_L.get_ibrav_type() == HEXAGONAL2) || 
                   (Rmg_L.get_ibrav_type() == ORTHORHOMBIC_PRIMITIVE) || 
                   (Rmg_L.get_ibrav_type() == CUBIC_PRIMITIVE) ||
                   (Rmg_L.get_ibrav_type() == MONOCLINIC_PRIMITIVE) ||
                   (Rmg_L.get_ibrav_type() == TRICLINIC_PRIMITIVE) ||
                   (Rmg_L.get_ibrav_type() == TETRAGONAL_PRIMITIVE));



    if(!special || (Rmg_L.get_ibrav_type() == HEXAGONAL) || 
                   (Rmg_L.get_ibrav_type() == HEXAGONAL2) ||
                   (Rmg_L.get_ibrav_type() == TRICLINIC_PRIMITIVE) ||
                   (Rmg_L.get_ibrav_type() == MONOCLINIC_PRIMITIVE))
        Rmg_T->trade_imagesx (a, rptr, dimx, dimy, dimz, images, FULL_TRADE);
    else
       Rmg_T->trade_imagesx (a, rptr, dimx, dimy, dimz, images, CENTRAL_TRADE);


    // Handle special combined operator first
    if(special && (order == APP_CI_EIGHT))
    {
        RmgTimer *RTA=NULL;
#if HIP_ENABLED || CUDA_ENABLED
        if(ct.use_gpu_fd)
        {
            if(ct.verbose) RTA = new RmgTimer("GPUFD");
            cc = FD.app8_combined(rptr, (DataType *)b, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec, true);
            if(ct.verbose) delete RTA;
        }
        else
        {
            if(ct.verbose) RTA = new RmgTimer("CPUFD");
            cc = FD.app8_combined (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec);
            if(ct.verbose) delete RTA;
        }
#else
        if(ct.verbose) RTA = new RmgTimer("CPUFD");
        cc = FD.app8_combined (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec);
        if(ct.verbose) delete RTA;
#endif

        return cc;
    }

    // First apply the laplacian
    if(order == APP_CI_SIXTH) {
        if(!special)
            cc = FiniteDiffLap (rptr, b, dimx, dimy, dimz, LC);
        else
            cc = FD.app6_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_EIGHT) {
        if(!special)
            cc = FiniteDiffLap (rptr, b, dimx, dimy, dimz, LC);
        else
            cc = FD.app8_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_TEN) {
        if(!special)
            cc = FiniteDiffLap (rptr, b, dimx, dimy, dimz, LC);
        else
            cc = FD.app10_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else {

        rmg_error_handler (__FILE__, __LINE__, "APP_DEL2 order not programmed yet in app_del2_driver.\n");
        return 0;   // Just to keep the compiler from complaining

    }

    if(!ct.is_gamma)
    {
        int pbasis = dimx*dimy*dimz;
        DataType *gx = new DataType[pbasis];
        DataType *gy = new DataType[pbasis];
        DataType *gz = new DataType[pbasis];
        std::complex<double> *kdr = new std::complex<double>[pbasis];


        if(special)
            FD.app_gradient_eighth (rptr, gx, gy, gz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        else
            FiniteDiffGrad((DataType *)rptr, gx, gy, gz, dimx, dimy, dimz, LC);

        // if complex orbitals compute dot products as well
        std::complex<double> I_t(0.0, 1.0);
        for(int idx = 0;idx < pbasis;idx++)
        {
            kdr[idx] = -I_t * (kvec[0] * std::complex<double>(std::real(gx[idx]), std::imag(gx[idx])) +
                               kvec[1] * std::complex<double>(std::real(gy[idx]), std::imag(gy[idx])) +
                               kvec[2] * std::complex<double>(std::real(gz[idx]), std::imag(gz[idx])));
        }
        if(typeid(DataType) == typeid(std::complex<double>))
        {
            std::complex<double> *pptr = (std::complex<double> *)b;
            for(int idx = 0;idx < pbasis;idx++) pptr[idx] = pptr[idx] - 2.0*kdr[idx];
        }
        if(typeid(DataType) == typeid(std::complex<float>))
        {
            std::complex<float> *pptr = (std::complex<float> *)b;
            for(int idx = 0;idx < pbasis;idx++) pptr[idx] = pptr[idx] - 2.0*kdr[idx];
        }
        
        delete [] kdr;
        delete [] gz;
        delete [] gy;
        delete [] gx;
    }

    return cc;

}


// The following two functions are for gamma point only
template <typename DataType>
double ApplyAOperator (DataType *a, DataType *b)
{
    int density = 1;

    int dimx = Rmg_G->get_PX0_GRID(density);
    int dimy = Rmg_G->get_PY0_GRID(density);
    int dimz = Rmg_G->get_PZ0_GRID(density);

    double gridhx = Rmg_G->get_hxgrid(density);
    double gridhy = Rmg_G->get_hygrid(density);
    double gridhz = Rmg_G->get_hzgrid(density);
                                                              
    return ApplyAOperator (&Rmg_L, Rmg_T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order);

}

template <typename DataType>
double ApplyAOperator (Lattice *L, TradeImages *T, DataType *a, DataType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    if(ct.kohn_sham_ke_fft)
    {
        FftLaplacianCoarse(a, b);
        FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
        DataType *ptr = NULL;
        // When ptr=NULL this does not do the finite differencing but just
        // returns the value of the diagonal element.
        double fd_diag = FD.app8_del2 (ptr, ptr, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        return fd_diag;
    }

    double cc = CPP_app_del2_driver (L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order, ct.alt_laplacian);
    return cc;
    
}
