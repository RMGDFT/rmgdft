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


/*

  This double Vdw::calculates the non-local correlation contribution to the
  energy and potential according to
 
     M. Dion, H. Rydberg, E. Schroeder, D. C. Langreth, and
     B. I. Lundqvist, Phys. Rev. Lett. 92, 246401 (2004).
 
  henceforth referred to as DION. Further information about the
  functional and its corresponding potential can be found in:
 
     T. Thonhauser, V.R. Cooper, S. Li, A. Puzder, P. Hyldgaard,
     and D.C. Langreth, Phys. Rev. B 76, 125112 (2007).
 
  The proper spin extension of vdW-DF, i.e. svdW-DF, is derived in
 
     T. Thonhauser, S. Zuluaga, C.A. Arter, K. Berland, E. Schroder,
     and P. Hyldgaard, Phys. Rev. Lett. 115, 136402 (2015).
 
  henceforth referred to as THONHAUSER.
 
 
  Two review article show many of the vdW-DF applications:
 
     D. C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009).
 
     K. Berland et al, Rep. Prog. Phys. 78, 066501 (2015).
 
 
  The method implemented is based on the method of G. Roman-Perez and
  J. M. Soler described in:
 
     G. Roman-Perez and J. M. Soler, PRL 103, 096102 (2009).
 
  henceforth referred to as SOLER.
 
  This code is based off of the vdW routines from Quantum Espresso.

*/ 

// Parallel fft library required for now

#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "vdW.h"
#include "RmgException.h"
#include "RmgSumAll.h"
#include "transition.h"
#include "packfuncs.h"
#include "RmgParallelFft.h"
#include "fft3d.h"
#if USE_LIBXC
#include "xc.h"
#endif

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
int Vdw::Nqs = VDW_NQPOINTS;
int Vdw::Nrpoints = VDW_NRPOINTS;
double Vdw::r_max;
bool Vdw::initialized = false;
double Vdw::q_mesh[VDW_NQPOINTS] = {
1.0e-5             , 0.0449420825586261, 0.0975593700991365, 0.159162633466142,
0.231286496836006, 0.315727667369529 , 0.414589693721418 , 0.530335368404141,
0.665848079422965, 0.824503639537924 , 1.010254382520950 , 1.227727621364570,
1.482340921174910, 1.780437058359530 , 2.129442028133640 , 2.538050036534580,
3.016440085356680, 3.576529545442460 , 4.232271035198720 , 5.0 };

double Vdw::kernel[VDW_NRPOINTS+1][VDW_NQPOINTS][VDW_NQPOINTS];
double Vdw::d2phi_dk2[VDW_NRPOINTS+1][VDW_NQPOINTS][VDW_NQPOINTS];

const double Vdw::epsr = 1.0e-12;
double Vdw::gmax;
double Vdw::dk;
double *Vdw::d2y_dx2;

/*

  rho_valence - valence charge
  rho_core    - core correction if applicable
  

*/
Vdw::Vdw (BaseGrid &G, Lattice &L, TradeImages &T, int type, double *rho_valence, double *rho_core, double &etxc, double &vtxc, double *v, bool gamma_flag)
{
  bool use_coarsegrid = !ct.use_vdwdf_finegrid;

  if((type != 1) && (type != 2))
      throw RmgFatalException() << "Vdw type was " << type << " but only 1 or 2 is allowed. " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

  // Grid parameters
  this->type = type;
  this->is_gamma = gamma_flag;

  // Gamma disabled until lookup table into fft grid is completed
  this->is_gamma = false;
  this->Grid = &G;
  this->T = &T;
  this->L = &L;
  this->pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
  this->hxgrid = G.get_hxgrid(G.default_FG_RATIO);
  this->hygrid = G.get_hygrid(G.default_FG_RATIO);
  this->hzgrid = G.get_hzgrid(G.default_FG_RATIO);
  this->dimx = G.get_PX0_GRID(G.default_FG_RATIO);
  this->dimy = G.get_PY0_GRID(G.default_FG_RATIO);
  this->dimz = G.get_PZ0_GRID(G.default_FG_RATIO);
  this->densgrid[0] = G.get_NX_GRID(G.default_FG_RATIO);
  this->densgrid[1] = G.get_NY_GRID(G.default_FG_RATIO);
  this->densgrid[2] = G.get_NZ_GRID(G.default_FG_RATIO);
  this->N = this->densgrid[0] * this->densgrid[1] * this->densgrid[2];

  this->pbasis_c = G.get_P0_BASIS(1);
  this->coarsegrid[0] = G.get_NX_GRID(1);
  this->coarsegrid[1] = G.get_NY_GRID(1);
  this->coarsegrid[2] = G.get_NZ_GRID(1);
  this->N_c = this->coarsegrid[0] * this->coarsegrid[1] * this->coarsegrid[2];

  // How many terms to include in the sum of SOLER equation 5.
  this->m_cut = 12;

  // Largest value of q_cut
  this->q_cut = q_mesh[Nqs-1];

  // smallest value of q_min
  this->q_min = q_mesh[0];

  // Local storage
  double *total_rho = new double[this->pbasis];
  double *gx = new double[3*this->pbasis];
  double *gy = gx + this->pbasis;
  double *gz = gy + this->pbasis;
  double *q0 = new double[this->pbasis]();
  double *dq0_drho = new double[this->pbasis]();
  double *dq0_dgradrho = new double[this->pbasis]();
  std::complex<double> *thetas = new std::complex<double> [this->pbasis*Nqs]();


  // Set up stuff that determines the precision of the intermediate calculations
  // if the use_coarsegrid is false then these are just the dimensions of the
  // input density grid. If use_coarsegrid is true then the dimensions of the
  // wavefunction grid are used.
  int calc_basis = this->pbasis;
  int N_calc = this->N;
  double *calc_rho=NULL, *calc_gx=NULL, *calc_gy=NULL, *calc_gz=NULL;
  Pw *planewaves_calc;
  if(use_coarsegrid) {
      calc_basis = this->pbasis_c;    // basis size for rho on this PE
      N_calc = this->N_c;             // total basis size across all nodes
      calc_rho = new double[calc_basis];
      calc_gx = new double[3*calc_basis];
      calc_gy = calc_gx + calc_basis;
      calc_gz = calc_gy + calc_basis;
  }

  // FFT plans
  if(use_coarsegrid) {

      this->plan_forward = fft_forward_coarse;
      this->plan_back = fft_backward_coarse;

  }
  else {

      this->plan_forward = fft_forward_fine;
      this->plan_back = fft_backward_fine;

  }

  // If not initialized we must read in the kernel table
  if(!this->initialized) {

      std::ifstream kernel_file;
      kernel_file.open("vdW_kernel_table", std::ios_base::in);
      if (!kernel_file)  { 
          throw RmgFatalException() << "Unable to open vdW_kernel_table file. Please make sure this file is in the current directory. " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
      }
     
      // Read in Nqs and Nrpoints from header and make sure it matches
      // these are the number of q points and the number of r points used for this kernel file, 
      std::string line;
      std::getline(kernel_file, line);
      std::istringstream is0(line);
      is0 >> Nqs >> Nrpoints;

      if((Nqs != VDW_NQPOINTS) || (Nrpoints != VDW_NRPOINTS)) {
          throw RmgFatalException() << "Mismatch for Nqs or Nrpoints, vdW_kernel_table appears to be corrupted." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
      }

      // Read in r_max the maximum allowed value of an r point
      std::getline(kernel_file, line);
      std::istringstream is1(line);
      is1 >> Vdw::r_max;
      Vdw::dk = 2.0 * PI / Vdw::r_max;
      
      // Read in the values of the q points used to generate this kernel.
      int idx = 0;
      for(int meshlines = 0;meshlines < 5;meshlines++) {
          std::getline(kernel_file, line);
          std::istringstream is2(line);
       
          while(is2 >> q_mesh[idx]) {
              idx++;
              if(idx == Nqs) break;
          }
      }

      // For each pair of q values, read in the function phi_q1_q2(k).
      // That is, the fourier transformed kernel function assuming q1 and
      // q2 for all the values of r used.
      int tlines = (Nrpoints + 1) / 4;   // 4 per line
      if((Nrpoints + 1) % 4) tlines++;
      
      for(int q1_i = 0;q1_i < Nqs;q1_i++) {
          for(int q2_i = 0;q2_i <= q1_i;q2_i++) {

              int idx = 0;
              for(int lines1 = 0;lines1 < tlines;lines1++) {
                  std::getline(kernel_file, line);
                  std::istringstream is3(line);

                  while(is3 >> kernel[idx][q1_i][q2_i]) {
                      idx++;
                      if(idx == (Nrpoints + 1)) break;
                  }
              }

              for(int jdx = 0;jdx < (Nrpoints + 1);jdx++) {
                  kernel[jdx][q2_i][q1_i] = kernel[jdx][q1_i][q2_i];
              }
          }
      }

      // Again, for each pair of q values (q1 and q2), read in the value of
      // the second derivative of the above mentiond Fourier transformed
      // kernel function phi_alpha_beta(k). These are used for spline
      // interpolation of the Fourier transformed kernel.
      for(int q1_i = 0;q1_i < Nqs;q1_i++) {
          for(int q2_i = 0;q2_i <= q1_i;q2_i++) {

              int idx = 0;
              for(int lines1 = 0;lines1 < tlines;lines1++) {
                  std::getline(kernel_file, line);
                  std::istringstream is4(line);

                  while(is4 >> d2phi_dk2[idx][q1_i][q2_i]) {
                      idx++;
                      if(idx == (Nrpoints + 1)) break;
                  }
              }

              for(int jdx = 0;jdx < (Nrpoints + 1);jdx++) {
                  d2phi_dk2[jdx][q2_i][q1_i] = d2phi_dk2[jdx][q1_i][q2_i];
              }

          }
      }

      kernel_file.close();

      // Allocate memory for the second derivatives used in the spline interpolation and initialize the values
      this->d2y_dx2 = new double[Nqs*Nqs]();
      double *y = new double[Nqs];
      double *temp_array = new double[Nqs];

      for(int P_i=0;P_i < Nqs;P_i++) {

          for(int ix=0;ix < Nqs;ix++) y[ix] = 0.0;          
          y[P_i] = 1.0;
          d2y_dx2[P_i] = 0.0;
          temp_array[0] = 0.0;

          for(int idx = 1;idx < Nqs - 1;idx++) {
              double temp1 = (q_mesh[idx] - q_mesh[idx-1]) / (q_mesh[idx+1] - q_mesh[idx-1]);
              double temp2 = temp1 * d2y_dx2[P_i + (idx-1)*Nqs] + 2.0;
              d2y_dx2[P_i + idx*Nqs] = (temp1 - 1.0) / temp2;

              temp_array[idx] = (y[idx+1] - y[idx]) / (q_mesh[idx+1] - q_mesh[idx]) - 
                                (y[idx] - y[idx-1]) / (q_mesh[idx] - q_mesh[idx-1]);
              temp_array[idx] = (6.0*temp_array[idx] / (q_mesh[idx+1] - q_mesh[idx-1]) - 
                                 temp1*temp_array[idx-1]) / temp2;
              
          }

          d2y_dx2[P_i + Nqs*(Nqs-1)] = 0.0;

          for(int idx=Nqs-2;idx >= 0;idx--) {
              d2y_dx2[P_i + idx*Nqs] = d2y_dx2[P_i + idx*Nqs] * d2y_dx2[P_i + (idx+1)*Nqs] + temp_array[idx];
          }

      }

      delete [] temp_array;
      delete [] y;

      this->initialized = true;
      
      // Print out run parameters and references
      this->info();
  }

  // Use the global plane wave objects
  if(use_coarsegrid) {
      planewaves_calc = coarse_pwaves;
  }
  else {
      planewaves_calc = fine_pwaves;
  }

  // Get total charge and compute it's gradient
  for(int i = 0;i < this->pbasis;i++) total_rho[i] = rho_valence[i] + rho_core[i];

  CPP_app_grad_driver (&L, &T, total_rho, gx, gy, gz, this->dimx, this->dimy, this->dimz, this->hxgrid, this->hygrid, this->hzgrid, APP_CI_TEN);

  // Have to generate half density versions of gradient and rho if use_coarsegrid is true.
  if(use_coarsegrid) {
      GetVtotPsi (calc_rho, total_rho, G.default_FG_RATIO);
      GetVtotPsi (calc_gx, gx, G.default_FG_RATIO);
      GetVtotPsi (calc_gy, gy, G.default_FG_RATIO);
      GetVtotPsi (calc_gz, gz, G.default_FG_RATIO);
  }
  else {
      calc_rho = total_rho;
      calc_gx = gx;
      calc_gy = gy;
      calc_gz = gz;
  }



  // --------------------------------------------------------------------
  // Find the value of q0 for all assigned grid points. q is defined in
  // equations 11 and 12 of DION and q0 is the saturated version of q
  // defined in equation 5 of SOLER. This routine also returns the
  // derivatives of the q0s with respect to the charge-density and the
  // gradient of the charge-density. These are needed for the potential
  // calculated below. This routine also calculates the thetas.

  this->get_q0_on_grid (calc_rho, q0, dq0_drho, dq0_dgradrho, thetas, calc_basis, calc_gx, calc_gy, calc_gz);

  double Ec_nl = this->vdW_energy(q0, thetas, calc_basis, N_calc, planewaves_calc);
  etxc += Ec_nl;



  // --------------------------------------------------------------------
  // Here we calculate the potential. This is calculated via equation 10
  // of SOLER, using the u_i(r) calculated from quations 11 and 12 of
  // SOLER. Each processor allocates the array to be the size of the full
  // grid  because, as can be seen in SOLER equation 10, processors need
  // to access grid points outside their allocated regions. Begin by
  // FFTing the u_i(k) to get the u_i(r) of SOLER equation 11.


  for(int iq = 0;iq < Nqs;iq++) {
    //pfft_execute_dft(plan_back_calc, (double (*)[2])&thetas[iq*calc_basis], (double (*)[2])&thetas[iq*calc_basis]);
      fft_3d((FFT_DATA *)&thetas[iq*calc_basis], (FFT_DATA *)&thetas[iq*calc_basis], 1, plan_back);
  }

  double *potential = new double[this->pbasis]();
  double *calc_potential = new double[this->pbasis]();
  this->get_potential(q0, dq0_drho, dq0_dgradrho, calc_potential, thetas, calc_basis, N_calc, calc_gx, calc_gy, calc_gz, planewaves_calc);

  if(use_coarsegrid) {
      FftInterpolation (G, calc_potential, potential, G.default_FG_RATIO);
  }
  else {
      for(int ix = 0;ix < this->pbasis;ix++)potential[ix] = calc_potential[ix];
  }

  // --------------------------------------------------------------------
  // The integral of rho(r)*potential(r) for the vtxc output variable.
  vtxc = 0.0;
  for(int ix=0;ix<this->pbasis;ix++) vtxc += rho_valence[ix] * potential[ix];
//  vtxc = RmgSumAll(vtxc, this->T->get_MPI_comm());
  vtxc = vtxc * L.omega / (double)this->N;
  
  for(int ix = 0;ix < this->pbasis;ix++) v[ix] += potential[ix];

  if(use_coarsegrid) {
      delete [] calc_gx;
      delete [] calc_rho;
  }

  delete [] calc_potential;
  delete [] potential;
  delete [] thetas;
  delete [] dq0_dgradrho;
  delete [] dq0_drho;
  delete [] q0;
  delete [] gx;
  delete [] total_rho;
}


Vdw::~Vdw(void)
{


}


// ####################################################################
//                         |                  |
//                         |  GET_Q0_ON_GRID  |
//                         |__________________|
//  
//  This routine first calculates the q value defined in (DION equations
//  11 and 12), then saturates it according to (SOLER equation 5). More
//  specifically it calculates the following:
//  
//       q0(ir) = q0 as defined above
//       dq0_drho(ir) = total_rho * d q0 /d rho
//       dq0_dgradrho = total_rho / |grad_rho| * d q0 / d |grad_rho|

void Vdw::get_q0_on_grid (double *calc_rho, double *q0, double *dq0_drho, double *dq0_dgradrho, std::complex<double> * thetas, int ibasis,
                          double *gx, double *gy, double *gz)
{


  double ec, dq0_dq;

  // Initialize q0-related arrays.
  for(int ix = 0;ix < ibasis;ix++) {
      q0[ix] = this->q_cut;      
      dq0_drho[ix] = 0.0;
      dq0_dgradrho[ix] = 0.0;
  }

  for(int ix = 0;ix < ibasis;ix++) {

      double trho = calc_rho[ix];

      // -----------------------------------------------------------------
      // This prevents numerical problems. If the charge density is
      // negative (an unphysical situation), we simply treat it as very
      // small. In that case, q0 will be very large and will be saturated.
      // For a saturated q0 the derivative dq0_dq will be 0 so we set q0 =
      // q_cut and dq0_drho = dq0_dgradrho = 0 and go on to the next
      // point.
      if(trho < this->epsr) continue; 

      // -----------------------------------------------------------------
      // Calculate some intermediate values needed to find q.

      double r_s = pow( 3.0 / (4.0*PI*trho) , (1.0/3.0));

      double s = sqrt( gx[ix]*gx[ix] + gy[ix]*gy[ix] + gz[ix]*gz[ix] ) /
         (2.0 * kF(trho) * trho );

      // -----------------------------------------------------------------
      // This is the q value defined in equations 11 and 12 of DION.
      // Use pw() from flib/functionals.f90 to get qc = kf/eps_x * eps_c.

      pw(r_s, 1, ec, dq0_drho[ix]);
      double q = -4.0*PI/3.0 * ec + kF(trho) * Fs(s);


      // -----------------------------------------------------------------
      // Bring q into its proper bounds.

       saturate_q ( q, q_cut, q0[ix], dq0_dq );
       if (q0[ix] < this->q_min) q0[ix] = this->q_min;


      // -----------------------------------------------------------------
      // Here we find derivatives. These are actually the density times
      // the derivative of q0 with respect to rho and grad_rho. The
      // density factor comes in since we are really differentiating
      // theta = (rho)*P(q0) with respect to density (or its gradient)
      // which will be
      //
      //    dtheta_drho = P(q0) + dP_dq0 * [rho * dq0_dq * dq_drho]
      //
      // and
      //
      //    dtheta_dgrad_rho = dP_dq0 * [rho * dq0_dq * dq_dgrad_rho]
      //
      // The parts in square brackets are what is calculated here. The
      // dP_dq0 term will be interpolated later.

      dq0_drho[ix]     = dq0_dq * trho * ( -4.0*PI/3.0 * 
                            (dq0_drho[ix] - ec)/trho + dqx_drho(trho, s) );
      dq0_dgradrho[ix] = dq0_dq * trho * kF(trho) * dFs_ds(s) * ds_dgradrho(trho);

  }

  // --------------------------------------------------------------------
  // Here we calculate the theta functions of SOLER equation 8. These are
  // defined as
  //
  //    rho * P_i(q0(rho, grad_rho))
  //
  // where P_i is a polynomial that interpolates a Kroneker delta
  // function at the point q_i (taken from the q_mesh) and q0 is the
  // saturated version of q. q is defined in equations 11 and 12 of DION
  // and the saturation proceedure is defined in equation 5 of SOLER.
  // This is the biggest memory consumer in the method since the thetas
  // array is (total # of FFT points)*Nqs complex numbers. In a parallel
  // run, each processor will hold the values of all the theta functions
  // on just the points assigned to it. thetas are stored in reciprocal
  // space as theta_i(k) because this is the way they are used later for
  // the convolution (equation 8 of SOLER). Start by interpolating the
  // P_i polynomials defined in equation 3 in SOLER for the particular q0
  // values we have.

  spline_interpolation (q_mesh, &Nqs, q0, &ibasis, thetas, d2y_dx2);

  for(int iq = 0;iq < Nqs;iq++) {
      for(int ix = 0;ix < ibasis;ix++) {
          thetas[ix + iq*ibasis] = thetas[ix + iq*ibasis] * calc_rho[ix];
      }
  }


  for(int iq = 0;iq < Nqs;iq++) {
      //pfft_execute_dft(plan_forward_calc, (double (*)[2])&thetas[iq*ibasis], (double (*)[2])&thetas[iq*ibasis]);
      fft_3d((FFT_DATA *)&thetas[iq*ibasis], (FFT_DATA *)&thetas[iq*ibasis], -1, plan_forward);
  }

}


double Vdw::vdW_energy(double *q0, std::complex<double> *thetas, int ibasis, int N_calc, Pw *pwaves)
{
  double *kernel_of_k = new double[Nqs*Nqs]();
  std::complex<double> *u_vdW = new std::complex<double>[Nqs * ibasis]();
  std::complex<double> *theta = new std::complex<double>[Nqs];

  double vdW_xc_energy = 0.0;
  double tpiba = 2.0 * PI / this->L->celldm[0];
  double G_multiplier = 1.0;
  if(is_gamma) G_multiplier = 2.0;

  for(int ig=0;ig < ibasis;ig++) {


      if(pwaves->gmask[ig] == 1.0){

          double g = sqrt(pwaves->gmags[ig]) * tpiba;

          this->interpolate_kernel(g, kernel_of_k);
          for(int idx=0;idx < Nqs;idx++) {
             theta[idx] = thetas[ig + idx*ibasis];
          }

          for(int q2_i=0;q2_i < Nqs;q2_i++) {

              for(int q1_i=0;q1_i < Nqs;q1_i++) {

                  u_vdW[q2_i*ibasis + ig] += kernel_of_k[q1_i*Nqs + q2_i] * theta[q1_i];

              }

              vdW_xc_energy = vdW_xc_energy + G_multiplier * std::real(u_vdW[q2_i*ibasis + ig] *
                              std::conj(theta[q2_i]));

          }
          // Special case for |g|=0 which is always the first g-vector on the
          // first node
          if((ig == 0) && (pct.gridpe == 0)) vdW_xc_energy /= G_multiplier;
      }
  }

  // Sum over all procs and normalize correctly. The fft of theta introduced a factor
  // of N which is squared by the product so we cancel it with 1/N^2. The volume element
  // for the reciprocal integral is (2pi)^3/omega which cancels one of the omega^2 terms
  // from Eq. 9 of Soler while the (2pi)^3 is cancelled by the implicit factor present
  // in k.
  // We skip the sum over processors for the returned quantity since that is the convention
  // followed in higher level routines where vdW_xc_energy is added to the other components
  // and then summed
  vdW_xc_energy = 0.5 * vdW_xc_energy * L->omega / (double)N_calc / (double)N_calc;
  double t1 = RmgSumAll(vdW_xc_energy, this->T->get_MPI_comm());
  rmg_printf("Van der Waals correlation energy = %16.9e Ha\n", t1);

  // Save u_vdW
  if(is_gamma) {
      // This is not correct. I need to conjugate the inverse vectors. Will need a lookup
      // table to map these correctly into the fft grid
      for(int ix=0;ix < ibasis*Nqs;ix++) thetas[ix] = std::conj(u_vdW[ix]);
  }
  else {
      for(int ix=0;ix < ibasis*Nqs;ix++) thetas[ix] = u_vdW[ix];
  }

  delete [] theta;
  delete [] u_vdW;
  delete [] kernel_of_k;
  return vdW_xc_energy;
}

void Vdw::get_potential(double *q0, double *dq0_drho, double *dq0_dgradrho, double *potential, std::complex<double> *u_vdW, 
                        int ibasis, int N_calc, double *gx, double *gy, double *gz, Pw *pwaves)
{

  std::complex<double> i(0.0,1.0);
  int q_low, q_hi, q;
  double *h_prefactor = new double[ibasis]();
  double *y = new double[Nqs]; 
  std::complex<double> *h = new std::complex<double>[ibasis]();
  double P, dP_dq0;
  double tpiba = 2.0 * PI / this->L->celldm[0];


  for(int ig = 0;ig < ibasis;ig++) {
     
      q_low = 0;
      q_hi = Nqs - 1;

      // -----------------------------------------------------------------
      // Figure out which bin our value of q0 is in in the q_mesh.
      while ( (q_hi - q_low) > 1) {

          q = ((q_hi + q_low)/2);

          if (q_mesh[q] > q0[ig]) {
             q_hi = q;
          }
          else {
             q_low = q;
          }

      } // end while

      if(q_low == q_hi)
           throw RmgFatalException() << "qhi == qlow" << __FILE__ << " at line " << __LINE__ << "\n";

      double dq = q_mesh[q_hi] - q_mesh[q_low];

      double a = (q_mesh[q_hi] - q0[ig]) / dq;
      double b = (q0[ig] - q_mesh[q_low]) / dq;
      double c = (a*a*a - a)*dq*dq/6.0;
      double d = (b*b*b - b)*dq*dq/6.0;
      double e = (3.0*a*a - 1.0)*dq/6.0;
      double f = (3.0*b*b - 1.0)*dq/6.0;

      for(int P_i = 0;P_i < Nqs;P_i++) {

          for(int iy = 0;iy < Nqs;iy++) y[iy] = 0.0;
          y[P_i] = 1.0;

          P      = a*y[q_low] + b*y[q_hi]  + c*d2y_dx2[P_i + q_low*Nqs] + d*d2y_dx2[P_i + q_hi*Nqs];
          dP_dq0 = (y[q_hi] - y[q_low])/dq - e*d2y_dx2[P_i + q_low*Nqs] + f*d2y_dx2[P_i + q_hi*Nqs];

         // --------------------------------------------------------------
         // The first term in equation 10 of SOLER.
         potential[ig] = potential[ig] + std::real(u_vdW[ig + P_i * ibasis] * (P + dP_dq0 * dq0_drho[ig]));
         if (q0[ig] != q_mesh[Nqs-1]) {
            h_prefactor[ig] = h_prefactor[ig] + std::real(u_vdW[ig + P_i * ibasis] * dP_dq0*dq0_dgradrho[ig]);
         }
      }
  }


  for(int icar = 0;icar < 3;icar++) {

     double *grad_rho;
     if(icar == 0) grad_rho = gx;
     if(icar == 1) grad_rho = gy;
     if(icar == 2) grad_rho = gz;

     for(int ix=0;ix < ibasis;ix++) {
         h[ix] = std::complex<double>(h_prefactor[ix] * grad_rho[ix], 0.0);
     }

     for(int ix = 0;ix < ibasis;ix++) {
        double gradient2 = gx[ix]*gx[ix] + gy[ix]*gy[ix] + gz[ix]*gz[ix];
        if ( gradient2 > 0.0) h[ix] = h[ix] / sqrt( gradient2 );
     }


     //pfft_execute_dft(plan_forward_calc, (double (*)[2])h, (double (*)[2])h);
     fft_3d((FFT_DATA *)h, (FFT_DATA *)h, -1, plan_forward);

     for(int ix=0;ix < ibasis;ix++) {
         if(pwaves->gmask[ix] == 1.0) {
             h[ix] = i * tpiba * pwaves->g[ix].a[icar] * h[ix];
         }
         else {
             h[ix] = std::complex<double>(0.0, 0.0);
         }
     }

     if (is_gamma) {
         // Not correct yet. Need a lookup table to map back into the fft grid
         for(int ix=0;ix < ibasis;ix++) {
             h[ix] = std::conj(h[ix]);
         }
     }

     //pfft_execute_dft(plan_back_calc, (double (*)[2])h, (double (*)[2])h);
     fft_3d((FFT_DATA *)h, (FFT_DATA *)h, 1, plan_back);
     double scale = 1.0 / (double)N_calc;

     for(int ix=0;ix < ibasis;ix++) potential[ix] -= scale * std::real(h[ix]);

  }


  // Now correct for the factor of N introduced by the earlier fft

  double scale = 1.0 / (double)N_calc;
  for(int ix=0;ix < ibasis;ix++) potential[ix] *= scale;

  //double echeck = 0.0;
  //for(int idx=0;idx<this->pbasis;idx++) echeck += this->total_rho[idx] * potential[idx];
  //echeck = RmgSumAll(echeck, this->T->get_MPI_comm());
  //rmg_printf("ECHECK = %18.8e\n",L->omega * echeck / (double)this->N);


  delete [] h;
  delete [] y;
  delete [] h_prefactor;
}
  
// ####################################################################
//                           |             |
//                           |  functions  |
//                           |_____________|
//
// Functions to be used in get_q0_on_grid and get_q0_on_grid_spin().

double Vdw::Fs(double s)
{
  double rFs, Z_ab;
  double fa=0.1234, fb=17.33, fc=0.163;  // Reparameterized values
                                         // from JCTC 5, 2745 (2009).

  if(this->type == 4) {
      rFs = pow( 1 + 15.0*fa*(s*s) + fb*pow(s,4) + fc*pow(s,6) ,(1.0/15.0));
  }
  else 
  {
      // ------------------------------------------------------------
      // Original functional choice for Fs, as definded in DION
      if (this->type == 1) Z_ab = -0.8491;
      if (this->type == 2) Z_ab = -1.887;
      rFs = 1.0 - Z_ab * (s*s) / 9.0;

  }

  return rFs;
}



double Vdw::dFs_ds(double s)
{
     double rdFs_ds, Z_ab;
     double fa=0.1234, fb=17.33, fc=0.163; // Reparameterized values
                                           // from JCTC 5, 2745 (2009).

     if(this->type == 4) {
         rdFs_ds = ( 30.0*fa*s + 4.0*fb*pow(s,3.0) + 6.0*fc*pow(s,5) ) 
                / ( 15.0*pow(( 1.0 + 15.0*fa*(s*s) + fb*pow(s,4) + fc*pow(s, 6) ),(14.0/15.0)) );
     }
     else 
     {
         // ------------------------------------------------------------
         // Original functional choice for Fs, as definded in DION
         if (this->type == 1) Z_ab = -0.8491;
         if (this->type == 2) Z_ab = -1.887;
         rdFs_ds =  -2.0 * s * Z_ab / 9.0;
     }
     return rdFs_ds;
}


double Vdw::kF(double rho)
{
    return (double)pow(( 3.0 * PI*PI * rho ), (1.0/3.0));
}


double Vdw::dkF_drho(double rho)
{
    return (1.0/3.0) * kF(rho) / rho;
}


double Vdw::ds_drho(double rho, double s)
{
    return -s * ( dkF_drho(rho) / kF(rho) + 1.0 / rho );
}


double Vdw::ds_dgradrho(double rho)
{
    return 1.0 / (2.0 * kF(rho) * rho);
}


double Vdw::dqx_drho(double rho, double s)
{
    return dkF_drho(rho) * Fs(s) + kF(rho) * dFs_ds(s) * ds_drho(rho, s);
}

void Vdw::saturate_q(double q, double q_cut, double &q0, double &dq0_dq)
{


    // --------------------------------------------------------------------
    // Here, we calculate q0 by saturating q according to equation 5 of
    // SOLER. Also, we find the derivative dq0_dq needed for the
    // derivatives dq0_drho and dq0_dgradrh0 discussed below.

    double e      = 0.0;
    dq0_dq = 0.0;
    for(int ix = 1;ix <= this->m_cut;ix++) {
       dq0_dq = dq0_dq + pow((q/q_cut),ix-1);
       e = e + pow((q/q_cut), ix)/(double)ix;
    }
    q0     = q_cut*(1.0 - exp(-e));
    dq0_dq = dq0_dq * exp(-e);
}

//-----------------------------------------------------------------------
//     iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
//     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
//
void Vdw::pw(double rs, int iflag, double &ec, double &vc)
{

  double a = 0.031091;
  double b1 = 7.5957;
  double b2 = 3.5876;
  double c0 = a;
  double c1 = 0.046644;
  double c2 = 0.00664;
  double c3 = 0.01043;
  double d0 = 0.4335;
  double d1 = 1.4408;
  double a1[2]={0.21370, 0.026481};
  double b3[2]={1.6382, -0.46647};
  double b4[2] = {0.49294, 0.13354};
  double lnrs, rs12, rs32, rs2, om, dom, olog;
  // 
  // high- and low-density formulae implemented but not used in PW case
  // (reason: inconsistencies in PBE/PW91 functionals)
  // 
  if ((rs < 1.0) && (iflag == 2)) {
     // high density formula
     lnrs = log (rs);
     ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs;
     vc = c0 * lnrs - (c1 + c0 / 3.0) + 2.0 / 3.0 * c2 * rs * lnrs - (2.0 * c3 + c2) / 3.0 * rs;
  }
  else if ((rs > 100.0) && (iflag==2)) {
     // low density formula
     ec = - d0 / rs + d1 / pow(rs,1.5);
     vc = - 4.0 / 3.0 * d0 / rs + 1.5 * d1 / pow(rs,1.5);
  }
  else {
     iflag--;
     // interpolation formula
     rs12 = sqrt (rs);
     rs32 = rs * rs12;
     rs2 = rs*rs;
     om = 2.0 * a * (b1 * rs12 + b2 * rs + b3[iflag] * rs32 + b4[iflag] * rs2);
     dom = 2.0 * a * (0.5 * b1 * rs12 + b2 * rs + 1.5 * b3[iflag] * rs32 + 2.0 * b4[iflag] * rs2);
     olog = log (1.0 + 1.0 / om);
     ec = - 2.0 * a * (1.0 + a1[iflag] * rs) * olog;
     vc = - 2.0 * a * (1.0 + 2.0 / 3.0 * a1[iflag] * rs) * olog - 2.0 / 3.0 * a * (1.0 + a1[iflag] * rs) * dom / (om * (om + 1.0) );
  }

}


// This routine is modeled after an algorithm from "Numerical Recipes in C" by
// Cambridge University Press, page 97.  Adapted for Fortran and the problem at
// hand.  This function is used to find the Phi_alpha_beta needed for equations
// 8 and 11 of SOLER.
//
// k = Input value, the magnitude of the g-vector for the current point.
//
// kernel_of_k = output array (allocated outside this routine)
// that holds the interpolated value of the kernel
// for each pair of q points (i.e. the phi_alpha_beta
// of the Soler method.
void Vdw::interpolate_kernel(double k, double *kernel_of_k)
{

  // --------------------------------------------------------------------
  // Check to make sure that the kernel table we have is capable of
  // dealing with this value of k. If k is larger than
  // Nr_points*2*pi/r_max then we can't perform the interpolation. In
  // that case, a kernel file should be generated with a larger number of
  // radial points.
  if ( k >= (double)Nrpoints*Vdw::dk ) {

      //std::cout << "GGGGGGG  " <<  k  << "  "  << (double)Nrpoints*Vdw::dk << std::endl;
      throw RmgFatalException() << "k value requested is out of range in " << __FILE__ << " at line " << __LINE__ << "\n";

  }

  for(int ix = 0;ix < Nqs*Nqs;ix++) kernel_of_k[ix] = 0.0;

  // This integer division figures out which bin k is in since the kernel
  // is set on a uniform grid.
  int k_i = (int)(k / Vdw::dk);

  double k_r, k_ii;
  k_r = modf(k/Vdw::dk, &k_ii);
  if(fabs(k_r) < this->epsr) {

      for(int q1_i = 0;q1_i < Nqs;q1_i++) {
          for(int q2_i = 0;q2_i <= q1_i;q2_i++) {

               kernel_of_k[q1_i + q2_i*Nqs] = kernel[k_i][q1_i][q2_i];
               kernel_of_k[q2_i + q1_i*Nqs] = kernel[k_i][q2_i][q1_i];

          }
      }
      return;
  } 
  
  double dk2 = dk * dk;
  double A = (Vdw::dk * ((double)k_i+1.0) - k) / Vdw::dk;
  double B = (k - Vdw::dk*(double)k_i) / Vdw::dk;
  double C = (A*A*A - A) * dk2 / 6.0;
  double D = (B*B*B - B) * dk2 / 6.0;

  for(int q1_i = 0;q1_i < Nqs;q1_i++) {
      for(int q2_i = 0;q2_i <= q1_i;q2_i++) {
      
          kernel_of_k[q1_i + q2_i*Nqs] = A*kernel[k_i][q1_i][ q2_i] + B*kernel[k_i+1][q1_i][q2_i]
             +(C*d2phi_dk2[k_i][q1_i][q2_i] + D*d2phi_dk2[k_i+1][q1_i][q2_i]);

          kernel_of_k[q2_i + q1_i*Nqs] = kernel_of_k[q1_i + q2_i*Nqs];

      }
  }

}


void Vdw::info(void) {
  // --------------------------------------------------------------------
  // Here we output some of the parameters being used in the run. This is
  // important because these parameters are read from the
  // vdW_kernel_table file. The user should ensure that these are the
  // parameters they were intending to use on each run.

  rmg_printf("\n     ************************************************************************\n");
  rmg_printf("     *                                                                      *\n");
  rmg_printf("     * You are using vdW-DF, which was implemented by the Thonhauser group. *\n");
  rmg_printf("     * Please cite the following two papers that made this development      *\n");
  rmg_printf("     * possible and the two reviews that describe the various versions:     *\n");
  rmg_printf("     *                                                                      *\n");
  rmg_printf("     *   T. Thonhauser et al., PRL 115, 136402 (2015).                      *\n");
  rmg_printf("     *   T. Thonhauser et al., PRB 76, 125112 (2007).                       *\n");
  rmg_printf("     *   K. Berland et al., Rep. Prog. Phys. 78, 066501 (2015).             *\n");
  rmg_printf("     *   D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009). *\n");
  rmg_printf("     *                                                                      *\n");
  rmg_printf("     *                                                                      *\n");
  rmg_printf("     * If you are calculating the stress with vdW-DF, please also cite:     *\n");
  rmg_printf("     *                                                                      *\n");
  rmg_printf("     *   R. Sabatini et al., J. Phys.: Condens. Matter 24, 424209 (2012).   *\n");
  rmg_printf("     *                                                                      *\n");
  rmg_printf("     ************************************************************************\n\n");
  
  rmg_printf("     Carrying out vdW-DF run using the following parameters:\n");
  rmg_printf("     Nqs    = %d  Npoints = %d  r_max = %12.8f\n", Nqs, Nrpoints, r_max );
  rmg_printf("     qmesh  = %12.8f  %12.8f  %12.8f  %12.8f\n",q_mesh[0],q_mesh[1],q_mesh[2],q_mesh[3]);
  for(int idx=4;idx < Nqs;idx+=4)
      rmg_printf("              %12.8f  %12.8f  %12.8f  %12.8f\n",q_mesh[idx],q_mesh[idx+1],q_mesh[idx+2],q_mesh[idx+3]);
  rmg_printf("\n\n");

}

