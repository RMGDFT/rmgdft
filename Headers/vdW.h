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


#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Pw.h"
#include "rmgtypedefs.h"
#include "rmg_mangling.h"

#include "vdW_params.h"

// Needed to deal with some issues when calling f90 module function from C++
#define     spline_interpolation                RMG_FC_MODULE(vdw_splines,spline_interpolation,mod_VDW_SPLINES, SPLINE_INTERPOLATION) 
#define     initialize_spline_interpolation     RMG_FC_MODULE(vdw_splines,initialize_spline_interpolation, mod_VDW_SPLINES, INITIALIZE_SPLINE_INTERPOLATION)

#define     generate_vdw_kernel       RMG_FC_MODULE(vdw_kernel, generate_vdw_kernel, mod_VDW_KERNEL, GENERATE_VDW_KERNEL)

extern "C" {void spline_interpolation (double *x, const int *Nx, double *evaluation_points, const int *Ngrid_points, std::complex<double>  *values, double *d2y_dx2);}
extern "C" {void initialize_spline_interpolation (double *x, const int *Nx, double *d2y_dx2);}
extern "C" {void generate_vdw_kernel(double *kernel, double *d2phi_dk2,
            double *q_mesh, MPI_Fint *intra_image_comm, int *mpime, int *nproc);}

class Vdw {

private:

    // BaseGrid class
    BaseGrid *Grid;

    // TradeImages object to use
    TradeImages *T;

    // Lattice object
    Lattice *L;

    // Pointer to second derivatives used in the spline interpolations
    static double *d2y_dx2;

    //
    int type;

    // Gamma flag
    bool is_gamma;

    // Real space basis on this node and grid parameters
    int pbasis;
    int pbasis_c;

    // Total number of grid points on the dense and coarse grids
    // and the number used for the intermediate calculations
    int N;
    int N_c;

    double hxgrid;
    double hygrid;
    double hzgrid;
    int dimx;
    int dimy;
    int dimz;
    ptrdiff_t densgrid[3];         // for passing to fft routines
    ptrdiff_t coarsegrid[3];         // for passing to fft routines

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
    static double_3d_array kernel;
    static double_3d_array d2phi_dk2;

    static double gmax;  // Maximum magnitude of g-vector

    // largest value of q_mesh
    double q_cut;

    // smallest value of q_mesh
    double q_min;

    // cutoff value for density
    static const double epsr;

    Pw *pwaves;

    double Fs(double s);
    double dFs_ds(double s);
    double kF(double rho);
    double dkF_drho(double rho);
    double ds_drho(double rho, double s);
    double ds_dgradrho(double rho);
    double dqx_drho(double rho, double s);
    void get_q0_on_grid (double *total_rho, double *q0, double *dq0_drho, double *dq0_dgradrho, std::complex<double> *thetas, 
                         int ibasis, double *gx, double *gy, double *gz);
    void saturate_q(double q, double q_cut, double &q0, double &dq0_dq);
    void pw(double rs, int iflag, double &ec, double &vc);
    void interpolate_kernel(double k, double *kernel_of_k);
    void interpolate_Dkernel_Dk (double k, double *dkernel_of_dk);
    void stress_vdW_DF_gradient (double *total_rho, double *grad_rho, double *q0, double *dq0_drho,
                                     double *dq0_dgradrho, std::complex<double> *thetas, double *sigma);
    void stress_vdW_DF_kernel (double *total_rho, double *q0, std::complex<double> *thetas, double *sigma);

public:
    Vdw (BaseGrid &G, Lattice &L, TradeImages &T, int type, double *rho_valence, double *rho_core, double &etxc, double &vtxc, double *v, bool gamma_flag);
    ~Vdw(void);

    double vdW_energy(double *q0, std::complex<double> *thetas, int ibasis, int N_calc);
    void get_potential(double *q0, double *dq0_drho, double *dq0_dgradrho, double *potential, std::complex<double> *u_vdW, 
                       int ibasis, int N_calc, double *gx, double *gy, double *gz);

    void stress_vdW_DF (double *rho_valence, double *rho_core, int nspin, double *sigma);
    void info(void);


};


#endif
