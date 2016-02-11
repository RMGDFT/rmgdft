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

#ifndef RMG_vdW_H
#define RMG_vdW_H 1

#if PFFT_LIBS

#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Pw.h"
#if USE_PFFT
    #include "pfft.h"
#endif

#define VDW_NQPOINTS  20
#define VDW_NRPOINTS  1024

extern "C" {void __vdw_splines_MOD_spline_interpolation (double *x, const int *Nx, double *evaluation_points, const int *Ngrid_points, std::complex<double>  *values, double *d2y_dx2);}
extern "C" {void __vdw_splines_MOD_initialize_spline_interpolation (double *x, const int *Nx, double *d2y_dx2);}

// Needed to deal with some issues when calling f90 module function from C++
#define  spline_interpolation  __vdw_splines_MOD_spline_interpolation
#define  initialize_spline_interpolation  __vdw_splines_MOD_initialize_spline_interpolation

#ifdef __cplusplus

class Vdw {

private:

    // BaseGrid class
    BaseGrid *Grid;

    // TradeImages object to use
    TradeImages *T;

    // Lattice object
    Lattice *L;

    // Pointer to plane wave object
    static Pw *plane_waves;

    // Pointer to second derivatives used in the spline interpolations
    static double *d2y_dx2;

    //
    int type;

    // Gamma flag
    bool is_gamma;

    // Real space basis on this node and grid parameters
    int pbasis;

    // Total number of grid points
    int N;

    double hxgrid;
    double hygrid;
    double hzgrid;
    int dimx;
    int dimy;
    int dimz;
    ptrdiff_t densgrid[3];         // for passing to fft routines
    double *total_rho;
    double *gx;
    double *gy;
    double *gz;
    double *q0;
    double *dq0_drho;
    double *dq0_dgradrho;

    // How many terms to include in the sum of SOLER equation 5.
    int m_cut;

    // ----------------------------------------------------------------------
    // The next 2 parameters define the q mesh to be used in the vdW_DF code.
    // These are perhaps the most important to have set correctly. Increasing
    // the number of q points will DRAMATICALLY increase the memory usage of
    // the vdW_DF code because the memory consumption depends quadratically
    // on the number of q points in the mesh. Increasing the number of q
    // points may increase accuracy of the vdW_DF code, although, in testing
    // it was found to have little effect. The largest value of the q mesh
    // is q_cut. All values of q0 (DION equation 11) larger than this value
    // during a run will be saturated to this value using equation 5 of
    // SOLER. In testing, increasing the value of q_cut was found to have
    // little impact on the results, though it is possible that in some
    // systems it may be more important. Always make sure that the variable
    // Nqs is consistent with the number of q points that are actually in the
    // variable q_mesh. Also, do not set any q value to 0. This will cause an
    // infinity in the Fourier transform.
    //
    //
    //                CHANGE THESE VALUES AT YOUR OWN RISK

    static bool initialized;
    static int Nqs;
    static int Nrpoints;
    static double r_max;
    static double dk;
    static double q_mesh[VDW_NQPOINTS];
    static double kernel[VDW_NRPOINTS+1][VDW_NQPOINTS][VDW_NQPOINTS];
    static double d2phi_dk2[VDW_NRPOINTS+1][VDW_NQPOINTS][VDW_NQPOINTS];
    static double gmax;  // Maximum magnitude of g-vector

    // largest value of q_mesh
    double q_cut;

    // smallest value of q_mesh
    double q_min;

    // cutoff value for density
    static const double epsr;

#if USE_PFFT
    pfft_plan plan_forward;
    pfft_plan plan_back;
#endif

public:
    Vdw (BaseGrid &G, Lattice &L, TradeImages &T, int type, double *rho_valence, double *rho_core, double &etxc, double &vtxc, double *v, bool gamma_flag);
    ~Vdw(void);

    void get_q0_on_grid (std::complex<double> *thetas);
    void saturate_q(double q, double q_cut, double &q0, double &dq0_dq);
    double vdW_energy(std::complex<double> *thetas);
    void get_potential(double *potential, std::complex<double> *u_vdW);

    double Fs(double s);
    double dFs_ds(double s);
    double kF(double rho);
    double dkF_drho(double rho);
    double ds_drho(double rho, double s);
    double ds_dgradrho(double rho);
    double dqx_drho(double rho, double s);
    void pw(double rs, int iflag, double &ec, double &vc);
    void index_to_gvector(int *index, double *gvec);
    void interpolate_kernel(double k, double *kernel_of_k);
    void fft_gradient(double *x, double *gx, double *gy, double *gz);
    void info(void);


};



#endif
#endif
#endif
