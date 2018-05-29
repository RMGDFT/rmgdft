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

// Applies Mehrstellen left hand side operator to a and returns result in b
// The first set of functions takes the input and output grids and a char string that defines
// the grid. More detailed parameters are then passed to the second set which may be accessed
// directly if more control is required.
//
// IN:    Input array a defined on coarse or fine grid
// OUT:   Output array b defined on coarse or fine grid
// IN:    grid = "Coarse" or "Fine" for grid type


template double ApplyAOperator<float>(float *, float *, float *, float *, float *, int, int, int, double, double, double, int);
template double ApplyAOperator<double>(double *, double *, double *, double *, double *,int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<float> >(std::complex<float> *, std::complex<float> *, std::complex<float> *, std::complex<float> *, std::complex<float> *,int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<double> >(std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *,int, int, int, double, double, double, int);

template double ApplyAOperator<float>(float *, float *, char *grid);
template double ApplyAOperator<double>(double *, double *, char *grid);
template double ApplyAOperator<std::complex<float> >(std::complex<float> *, std::complex<float> *, char *grid);
template double ApplyAOperator<std::complex<double> >(std::complex<double> *, std::complex<double> *, char *grid);

template double ApplyAOperator<float>(float *, float *, char *grid, BaseGrid *G, TradeImages *T);
template double ApplyAOperator<double>(double *, double *, char *grid, BaseGrid *G, TradeImages *T);
template double ApplyAOperator<std::complex<float> >(std::complex<float> *, std::complex<float> *, char *grid, BaseGrid *G, TradeImages *T);
template double ApplyAOperator<std::complex<double> >(std::complex<double> *, std::complex<double> *, char *grid, BaseGrid *G, TradeImages *T);

template double ApplyAOperator<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int);
template double ApplyAOperator<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);


template <typename DataType>
double ApplyAOperator (DataType *a, DataType *b, DataType *gx, DataType *gy, DataType *gz, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{
    
    double cc = 0.0;
    FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
    int sbasis = (dimx + order) * (dimy + order) * (dimz + order);
    int images = order / 2;
    size_t alloc = (sbasis + 64) * sizeof(DataType);
    DataType *rptr;

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


    Rmg_T->trade_imagesx (a, rptr, dimx, dimy, dimz, images, CENTRAL_TRADE);

    // First apply the laplacian
    if(order == APP_CI_FOURTH) {
        cc = FD.app4_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        if(!ct.is_gamma) FD.app_gradient_fourth (rptr, gx, gy, gz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_SIXTH) {
        cc = FD.app6_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        if(!ct.is_gamma) FD.app_gradient_sixth (rptr, gx, gy, gz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_EIGHT) {
        cc = FD.app8_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        if(!ct.is_gamma) FD.app_gradient_eighth (rptr, gx, gy, gz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_TEN) {
        cc = FD.app10_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        if(!ct.is_gamma) FD.app_gradient_tenth (rptr, gx, gy, gz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_TWELVE) {
        cc = FD.app12_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        if(!ct.is_gamma) FD.app_gradient_twelfth (rptr, gx, gy, gz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else {

        rmg_error_handler (__FILE__, __LINE__, "APP_DEL2 order not programmed yet in app_del2_driver.\n");
        return 0;   // Just to keep the compiler from complaining

    }

    if(alloc > 110592) delete [] rptr;

    return cc;

}

template <typename DataType>
double ApplyAOperator (DataType *a, DataType *b, char *grid)
{
    int density;
    const char *coarse = "Coarse";
    const char *fine = "Fine";

    if(!strcmp(grid, coarse)) {
        density = 1;
    }
    else if(!strcmp(grid, fine)) {
        density = Rmg_G->default_FG_RATIO;
    }
    else {
        throw RmgFatalException() << "Error! Grid type " << grid << " not defined in "
                                 << __FILE__ << " at line " << __LINE__ << "\n";
    }

    int dimx = Rmg_G->get_PX0_GRID(density);
    int dimy = Rmg_G->get_PY0_GRID(density);
    int dimz = Rmg_G->get_PZ0_GRID(density);

    double gridhx = Rmg_G->get_hxgrid(density);
    double gridhy = Rmg_G->get_hygrid(density);
    double gridhz = Rmg_G->get_hzgrid(density);
                                                              
    return ApplyAOperator (&Rmg_L, Rmg_T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order);

}

template <typename DataType>
double ApplyAOperator (DataType *a, DataType *b, char *grid, BaseGrid *G, TradeImages *T)
{
    int density;
    const char *coarse = "Coarse";
    const char *fine = "Fine";

    if(!strcmp(grid, coarse)) {
        density = 1;
    }
    else if(!strcmp(grid, fine)) {
        density = G->default_FG_RATIO;
    }
    else {
        throw RmgFatalException() << "Error! Grid type " << grid << " not defined in "
                                 << __FILE__ << " at line " << __LINE__ << "\n";
    }

    int dimx = G->get_PX0_GRID(density);
    int dimy = G->get_PY0_GRID(density);
    int dimz = G->get_PZ0_GRID(density);

    double gridhx = G->get_hxgrid(density);
    double gridhy = G->get_hygrid(density);
    double gridhz = G->get_hzgrid(density);
                                                              
    return ApplyAOperator (&Rmg_L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order);

}


template <typename DataType>
double ApplyAOperator (Lattice *L, TradeImages *T, DataType *a, DataType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    if(ct.discretization == MEHRSTELLEN_DISCRETIZATION) {

        return CPP_app_cil_driver (L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);

    }
    else if(ct.discretization == CENTRAL_DISCRETIZATION) {

        double cc = CPP_app_del2_driver (L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order, ct.alt_laplacian);
        return cc;

    }
    
    throw RmgFatalException() << "Error! Unknown discretization method " << " in "
                                 << __FILE__ << " at line " << __LINE__ << "\n";
    return 0.0;
}

