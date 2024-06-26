/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#ifndef RMG_FiniteDiff_H
#define RMG_FiniteDiff_H 1

/* Order of finite differencing for driver routines */
#define APP_CI_FFT 0
#define APP_CI_SECOND 2
#define APP_CI_FOURTH 4
#define APP_CI_SIXTH 6
#define APP_CI_EIGHT 8
#define APP_CI_TEN 10
#define APP_CI_TWELVE 12

#ifdef __cplusplus
#include <unordered_map>
#include "Lattice.h"
#include "TradeImages.h"
#include "boundary_conditions.h"
#include "GpuAlloc.h"
#include "LaplacianCoeff.h"

template <typename RmgType>
void CPP_app_cir_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order);
template <typename RmgType>
double CPP_app_cil_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);
template <typename RmgType>
void CPP_app_grad_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * bx, RmgType * by, RmgType * bz, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order, bool alt_flag);
template <typename RmgType>
double CPP_app_del2_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);
template <typename RmgType>
double CPP_app_del2_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order, bool alt_flag);
template <typename RmgType>
double CPP_app_del2_driver_int (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order, bool alt_flag);



#include "rmg_error.h"

class FiniteDiff {

private:
    Lattice *L;
    BaseGrid *G;

    int x_type;
    int y_type;
    int z_type;

    bool alt_laplacian;
    // For non-periodic boundary conditions.
    double *np_weight;
    double *np_xweight;
    double *np_yweight;
    double *np_zweight;
    int np_density;
    int *xoff;
    int *yoff;
    int *zoff;
    int order;
    int stride;


public:
    FiniteDiff(Lattice *lptr);
    FiniteDiff(Lattice *lptr, bool alt_flag);
    FiniteDiff(Lattice *lptr, BaseGrid *G, int xtype, int ytype, int ztype, int density, int order);
    static void gen_weights(int n, int m, double xr, double *x, double *w);
    static void set_allocation_limit(int lim);
    static int LCkey(double a0h);
    void set_alt_laplacian_flag(bool flag);
    static int allocation_limit;
    static double cfac[13];

    // Used to access Coeffs for a given grid and order.
    // The key is dimx*dimy*dimz+order
    static std::unordered_map<int, LaplacianCoeff *> FdCoeffs;


    ~FiniteDiff(void);

    bool check_anisotropy(double hx, double hy, double hz, double limit);

    template <typename RmgType>
    double app_del2_np (RmgType *rptr, RmgType *b, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app_cil_sixth (RmgType *rptr, RmgType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    void app_cir_sixth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz);

    template <typename RmgType>
    double app2_del2 (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app2_del2_offset (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int offset);

    template <typename RmgType>
    double app8_del2(RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app10_del2(RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app6_del2(RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app12_del2(RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app_cil_fourth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    void app_cir_fourth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz);

    template <typename RmgType>
    void app_gradient_sixth (RmgType * rptr, RmgType * wxr, RmgType *wyr, RmgType *wzr, int dimx, int dimy, int dimz,
                                   double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    void app_gradient_eighth (RmgType * rptr, RmgType * wxr, RmgType *wyr, RmgType *wzr, int dimx, int dimy, int dimz,
                                   double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    void app6_gradient_general (RmgType * __restrict__ a,
                                RmgType * __restrict__ gx,
                                RmgType * __restrict__ gy,
                                RmgType * __restrict__ gz,
                                int dimx, int dimy, int dimz);

    template <typename RmgType, int order>
    void fd_gradient_general (RmgType * __restrict__ a,
                                RmgType * __restrict__ gx,
                                RmgType * __restrict__ gy,
                                RmgType * __restrict__ gz,
                                double hxgrid,
                                int dimx, int dimy, int dimz);

    template <typename RmgType>
    void app_gradient_tenth (RmgType * rptr, RmgType * wxr, RmgType *wyr, RmgType *wzr, int dimx, int dimy, int dimz,
                                   double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app6_combined(
		    RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz,
		    double *kvec, bool use_gpu);

    template <typename RmgType, int order>
    double app_combined(
		    RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz,
		    double *kvec, bool use_gpu);

    double fd_coeff0(int order, double hxgrid);

    template <typename RmgType>
    void fd_combined_coeffs(int order, double hxgrid, int axis, RmgType * cm, RmgType *cp, double *kvec);

    template <typename RmgType>
    void fd_gradient_coeffs(int order, double hxgrid, int axis, RmgType *cx, RmgType *cy, RmgType *cz);

    double app6_coeff0(void);

    template <typename RmgType>
    void app6_combined_coeffs(int order, int ax, RmgType * cm, RmgType *cp, double *kvec);

    template <typename RmgType>
    void app6_gradient_coeffs(int order, int axis , RmgType *cx, RmgType *cy, RmgType *cz);

};
#endif

#endif
