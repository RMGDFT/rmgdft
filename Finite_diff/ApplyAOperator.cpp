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
#include "rmg_complex.h"

// Applies Mehrstellen left hand side operator to a and returns result in b
// The first set of functions takes the input and output grids and a char string that defines
// the grid. More detailed parameters are then passed to the second set which may be accessed
// directly if more control is required.
//
// IN:    Input array a defined on coarse or fine grid
// OUT:   Output array b defined on coarse or fine grid
// IN:    grid = "Coarse" or "Fine" for grid type


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
    double cc = 0.0;
    FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
    int sbasis = (dimx + order) * (dimy + order) * (dimz + order);
    int images = order / 2;
    size_t alloc = (sbasis + 64) * sizeof(DataType);
    DataType *rptr;
    int special = ((Rmg_L.get_ibrav_type() == HEXAGONAL) || (Rmg_L.get_ibrav_type() == ORTHORHOMBIC_PRIMITIVE) || (Rmg_L.get_ibrav_type() == CUBIC_PRIMITIVE));

    // while alloca is dangerous it's very fast for small arrays and the 110k limit
    // is fine for linux and 64bit power
    if(alloc <= 110592)
    {
        rptr = (DataType *)alloca(alloc);
    }
    else
    {
        rptr = new DataType[sbasis + 64];
    }



    if(!special || (Rmg_L.get_ibrav_type() == HEXAGONAL))
        Rmg_T->trade_imagesx (a, rptr, dimx, dimy, dimz, images, FULL_TRADE);
    else
       Rmg_T->trade_imagesx (a, rptr, dimx, dimy, dimz, images, CENTRAL_TRADE);


    // Handle special combined operator first
    if(special && (order == APP_CI_EIGHT) && (Rmg_L.get_ibrav_type() != HEXAGONAL))
    {
        cc = FD.app8_combined (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec);
        if(alloc > 110592) delete [] rptr;
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
            FiniteDiffGrad(rptr, gx, gy, gz, dimx, dimy, dimz, LC);

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

    if(alloc > 110592) delete [] rptr;

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
