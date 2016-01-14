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

#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>
#include <fstream>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "vdW.h"
#include "RmgException.h"
#include "xc.h"
#include "RmgSumAll.h"
#include "transition.h"
#include "Pw.h"
#include "pfft.h"

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

//const double Vdw::epsr = 1.0e-12;
const double Vdw::epsr = 1.0e-12;
double Vdw::gmax;
double Vdw::dk;
Pw *Vdw::plane_waves;
double *Vdw::d2y_dx2;

/*

  rho_valence - valence charge
  rho_core    - core correction if applicable
  

*/

Vdw::Vdw (BaseGrid &G, Lattice &L, TradeImages &T, int type, double *rho_valence, double *rho_core, double &etxc, double &vtxc, double *v)
{

  // Grid parameters
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

  // How many terms to include in the sum of SOLER equation 5.
  this->m_cut = 12;

  // Largest value of q_cut
  this->q_cut = Vdw::q_mesh[Vdw::Nqs-1];

  // smallest value of q_min
  this->q_min = Vdw::q_mesh[0];

  // Local storage
  this->total_rho = new double[this->pbasis];
  this->gx = new double[this->pbasis];
  this->gy = new double[this->pbasis];
  this->gz = new double[this->pbasis];
  this->q0 = new double[this->pbasis]();
  this->dq0_drho = new double[this->pbasis]();
  this->dq0_dgradrho = new double[this->pbasis]();
  this->thetas = new std::complex<double> [this->pbasis*Vdw::Nqs]();
  double *potential = new double[this->pbasis]();


  // If not initialized we must read in the kernel table
  if(!this->initialized) {

      std::ifstream kernel_file;
      kernel_file.open("vdW_kernel_table", std::ios_base::in);
      if (!kernel_file)  { 
          throw RmgFatalException() << "Unable to open vdW_kernel_table file. Please make sure this file is in the current directory. " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
      }
     
      // Read in Vdw::Nqs and Vdw::Nrpoints from header and make sure it matches
      // these are the number of q points and the number of r points used for this kernel file, 
      std::string line;
      std::getline(kernel_file, line);
      std::istringstream is0(line);
      is0 >> Vdw::Nqs >> Vdw::Nrpoints;

      if((Vdw::Nqs != VDW_NQPOINTS) || (Vdw::Nrpoints != VDW_NRPOINTS)) {
          throw RmgFatalException() << "Mismatch for Vdw::Nqs or Vdw::Nrpoints, vdW_kernel_table appears to be corrupted." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
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
       
          while(is2 >> Vdw::q_mesh[idx]) {
              idx++;
              if(idx == Vdw::Nqs) break;
          }
      }

      // For each pair of q values, read in the function phi_q1_q2(k).
      // That is, the fourier transformed kernel function assuming q1 and
      // q2 for all the values of r used.
      int tlines = (Vdw::Nrpoints + 1) / 4;   // 4 per line
      if((Vdw::Nrpoints + 1) % 4) tlines++;
      
      for(int q1_i = 0;q1_i < Vdw::Nqs;q1_i++) {
          for(int q2_i = 0;q2_i <= q1_i;q2_i++) {
              for(int lines1 = 0;lines1 < tlines;lines1++) {
                  std::getline(kernel_file, line);
                  int idx = 0;
                  std::istringstream is3(line);

                  while(is3 >> Vdw::kernel[idx][q1_i][q2_i]) {
                      idx++;
                      if(idx == (Vdw::Nrpoints + 1)) break;
                  }
              }

              for(int jdx = 0;jdx < (Vdw::Nrpoints + 1);jdx++) {
                  Vdw::kernel[jdx][q2_i][q1_i] = Vdw::kernel[jdx][q1_i][q2_i];
              }
          }
      }

      // Again, for each pair of q values (q1 and q2), read in the value of
      // the second derivative of the above mentiond Fourier transformed
      // kernel function phi_alpha_beta(k). These are used for spline
      // interpolation of the Fourier transformed kernel.
      for(int q1_i = 0;q1_i < Vdw::Nqs;q1_i++) {
          for(int q2_i = 0;q2_i <= q1_i;q2_i++) {

              int idx = 0;
              for(int lines1 = 0;lines1 < tlines;lines1++) {
                  std::getline(kernel_file, line);
                  std::istringstream is4(line);

                  while(is4 >> Vdw::d2phi_dk2[idx][q1_i][q2_i]) {
                      idx++;
                      if(idx == (Vdw::Nrpoints + 1)) break;
                  }
              }

              for(int jdx = 0;jdx < (Vdw::Nrpoints + 1);jdx++) {
                  Vdw::d2phi_dk2[jdx][q2_i][q1_i] = Vdw::d2phi_dk2[jdx][q1_i][q2_i];
              }

          }
      }

      kernel_file.close();

      // Plane wave object
      this->plane_waves = new Pw(G, L, G.default_FG_RATIO);

      // Allocate memory for the second derivatives used in the spline interpolation and initialize the values
      this->d2y_dx2 = new double[Vdw::Nqs*Vdw::Nqs];
      double *y = new double[Vdw::Nqs];
      double *temp_array = new double[Vdw::Nqs];

      for(int P_i=0;P_i < Vdw::Nqs;P_i++) {

          for(int ix=0;ix < Vdw::Nqs;ix++) y[ix] = 0.0;          
          y[P_i] = 1.0;
          d2y_dx2[P_i] = 0.0;
          temp_array[0] = 0.0;

          for(int idx = 1;idx < Vdw::Nqs - 1;idx++) {
              double temp1 = (Vdw::q_mesh[idx] - Vdw::q_mesh[idx-1]) / (Vdw::q_mesh[idx+1] - Vdw::q_mesh[idx-1]);
              double temp2 = temp1 * d2y_dx2[P_i + (idx-1)*Vdw::Nqs] + 2.0;
              d2y_dx2[P_i + idx*Vdw::Nqs] = (temp1 - 1.0) / temp2;

              temp_array[idx] = (y[idx+1] - y[idx]) / (Vdw::q_mesh[idx+1] - Vdw::q_mesh[idx]) - 
                                (y[idx] - y[idx-1]) / (Vdw::q_mesh[idx] - Vdw::q_mesh[idx-1]);
              temp_array[idx] = (6.0*temp_array[idx] / (Vdw::q_mesh[idx+1] - Vdw::q_mesh[idx-1]) - 
                                 temp1*temp_array[idx-1]) / temp2;
              
          }

          d2y_dx2[P_i + Vdw::Nqs*(Vdw::Nqs-1)] = 0.0;

          for(int idx=Vdw::Nqs-2;idx >= 1;idx--) {
              d2y_dx2[P_i + idx*Vdw::Nqs] = d2y_dx2[P_i + idx*Vdw::Nqs] * d2y_dx2[P_i + (idx+1)*Vdw::Nqs] + temp_array[idx];
          }

      }

      delete [] temp_array;
      delete [] y;

      this->initialized = true;
      
  }

  this->plan_forward = pfft_plan_dft_3d(this->densgrid, (double (*)[2])this->thetas, (double (*)[2])this->thetas,
                                                 pct.pfft_comm, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_DESTROY_INPUT | PFFT_ESTIMATE);

  
  // Get total charge and compute it's gradient
  for(int i = 0;i < pbasis;i++) total_rho[i] = rho_valence[i] + rho_core[i];

  CPP_app_grad_driver (&L, &T, total_rho, gx, gy, gz, this->dimx, this->dimy, this->dimz, this->hxgrid, this->hygrid, this->hzgrid, APP_CI_SIXTH);
  
  // --------------------------------------------------------------------
  // Find the value of q0 for all assigned grid points. q is defined in
  // equations 11 and 12 of DION and q0 is the saturated version of q
  // defined in equation 5 of SOLER. This routine also returns the
  // derivatives of the q0s with respect to the charge-density and the
  // gradient of the charge-density. These are needed for the potential
  // calculated below. This routine also calculates the thetas.

  this->get_q0_on_grid ();

  double Ec_nl = this->vdW_energy();
  etxc += Ec_nl;


  // Now get the potential
  this->plan_back = pfft_plan_dft_3d(this->densgrid, (double (*)[2])this->thetas, (double (*)[2])this->thetas,
                                                 pct.pfft_comm, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_DESTROY_INPUT | PFFT_ESTIMATE);
 
  this->get_potential(potential, this->thetas);
   
  for(int ix = 0;ix < this->pbasis;ix++) v[ix] += potential[ix];

  delete [] potential;
}


// Destructor just frees memory
Vdw::~Vdw(void)
{

  pfft_destroy_plan(this->plan_back);
  pfft_destroy_plan(this->plan_forward);
  delete [] this->thetas;
  delete [] this->dq0_dgradrho;
  delete [] this->dq0_drho;
  delete [] this->q0;
  delete [] this->gz;
  delete [] this->gy;
  delete [] this->gx;
  delete [] this->total_rho;

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

void Vdw::get_q0_on_grid (void)
{

  double ec, dq0_dq;

  // Initialize q0-related arrays.
  for(int ix = 0;ix < this->pbasis;ix++) {
      this->q0[ix] = this->q_cut;      
      this->dq0_drho[ix] = 0.0;
      this->dq0_dgradrho[ix] = 0.0;
  }

  for(int ix = 0;ix < this->pbasis;ix++) {

      double trho = this->total_rho[ix];

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

      double s = sqrt( this->gx[ix]*this->gx[ix] + this->gy[ix]*this->gy[ix] + this->gz[ix]*this->gz[ix] ) /
         (2.0 * kF(trho) * trho );

      // -----------------------------------------------------------------
      // This is the q value defined in equations 11 and 12 of DION.
      // Use pw() from flib/functionals.f90 to get qc = kf/eps_x * eps_c.

      pw(r_s, 1, ec, dq0_drho[ix]);
      double q = -4.0*PI/3.0 * ec + kF(trho) * Fs(s);


      // -----------------------------------------------------------------
      // Bring q into its proper bounds.

       saturate_q ( q, q_cut, q0[ix], dq0_dq );
       if (this->q0[ix] < this->q_min) this->q0[ix] = this->q_min;


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

      this->dq0_drho[ix]     = dq0_dq * trho * ( -4.0*PI/3.0 * 
                            (this->dq0_drho[ix] - ec)/trho + dqx_drho(trho, s) );
      this->dq0_dgradrho[ix] = dq0_dq * trho * kF(trho) * dFs_ds(s) * ds_dgradrho(trho);

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

  spline_interpolation (Vdw::q_mesh, &Vdw::Nqs, q0, &this->pbasis, thetas);

  for(int iq = 0;iq < Vdw::Nqs;iq++) {
      for(int ix = 0;ix < this->pbasis;ix++) {
          thetas[ix + iq*this->pbasis] = thetas[ix + iq*this->pbasis] * this->total_rho[ix];
      }
  }

  for(int iq = 0;iq < Vdw::Nqs;iq++) {
      pfft_execute_dft(plan_forward, (double (*)[2])&thetas[iq*this->pbasis], (double (*)[2])&thetas[iq*this->pbasis]);
  }
}


double Vdw::vdW_energy(void)
{
  double *kernel_of_k = new double[Vdw::Nqs*Vdw::Nqs]();
  std::complex<double> *u_vdW = new std::complex<double>[Vdw::Nqs * this->pbasis]();
  std::complex<double> *theta = new std::complex<double>[Vdw::Nqs];

  double vdW_xc_energy = 0.0;
  double tpiba = 2.0 * PI / this->L->celldm[0];
  for(int ig=0;ig < this->pbasis;ig++) {

      double g = this->plane_waves->gmags[ig] * tpiba;
      this->interpolate_kernel(g, kernel_of_k);
      for(int idx=0;idx < Vdw::Nqs;idx++) theta[idx] = thetas[ig + idx*Vdw::Nqs];

      for(int q2_i=0;q2_i < Vdw::Nqs;q2_i++) {

          for(int q1_i=0;q1_i < Vdw::Nqs;q1_i++) {

              u_vdW[q2_i*Vdw::Nqs + ig] += kernel_of_k[q1_i*Vdw::Nqs + q2_i] * theta[q1_i] * this->plane_waves->gmask[ig];

          }

          vdW_xc_energy = vdW_xc_energy + std::real(u_vdW[q2_i*Vdw::Nqs + ig] *
                          std::conj(theta[q2_i]));

      }

  }

  double t1 = RmgSumAll(vdW_xc_energy, this->T->get_MPI_comm());
  vdW_xc_energy = t1 * 0.5 * L->omega / (double)this->N;
  rmg_printf("Van der Waals correlation energy = %16.9f Ha\n", vdW_xc_energy);

  delete [] theta;
  delete [] u_vdW;
  delete [] kernel_of_k;
  return vdW_xc_energy;
}

void Vdw::get_potential(double *potential, std::complex<double> *u_vdW)
{

  std::complex<double> i(0.0,1.0);
  int q_low, q_hi, q;
  double *h_prefactor = new double[this->pbasis];
  double *y = new double[Vdw::Nqs]; 
  std::complex<double> *h = new std::complex<double>[this->pbasis];
  double P, dP_dq0;
  double tpiba = 2.0 * PI / this->L->celldm[0];

  for(int ig = 0;ig < this->pbasis;ig++) {
     
      q_low = 0;
      q_hi = Vdw::Nqs - 1; 

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

      double dq = Vdw::q_mesh[q_hi] - Vdw::q_mesh[q_low];

      double a = (Vdw::q_mesh[q_hi] - q0[ig]) / dq;
      double b = (q0[ig] - Vdw::q_mesh[q_low]) / dq;
      double c = (pow(a,3.0) - a)*pow(dq,2/6.0);
      double d = (pow(b,3.0) - b)*pow(dq,2/6.0);
      double e = (3.0*pow(a,2.0) - 1.0)*dq/6.0;
      double f = (3.0*pow(b,2.0) - 1.0)*dq/6.0;

      for(int P_i = 0;P_i < Vdw::Nqs;P_i++) {

          for(int iy = 0;iy < Vdw::Nqs;iy++) y[iy] = 0.0;
          y[P_i] = 1.0;

          P      = a*y[q_low] + b*y[q_hi]  + c*d2y_dx2[P_i + q_low*Vdw::Nqs] + d*Vdw::d2y_dx2[P_i + q_hi*Vdw::Nqs];
          dP_dq0 = y[q_hi] - y[q_low]/dq - e*d2y_dx2[P_i + q_low*Vdw::Nqs] + f*Vdw::d2y_dx2[P_i + q_hi*Vdw::Nqs];


         // --------------------------------------------------------------
         // The first term in equation 10 of SOLER.
         potential[ig] = potential[ig] + std::real(u_vdW[ig + P_i * this->pbasis] * (P + dP_dq0 * dq0_drho[ig]));
         if (q0[ig] != q_mesh[Vdw::Nqs-1]) {
            h_prefactor[ig] = h_prefactor[ig] + std::real(u_vdW[ig + P_i * this->pbasis] * dP_dq0*dq0_dgradrho[ig]);
         }
      }


  }


  for(int icar = 0;icar < 3;icar++) {

     double *grad_rho;
     if(icar == 0) grad_rho = this->gx;
     if(icar == 1) grad_rho = this->gy;
     if(icar == 2) grad_rho = this->gz;

     for(int ix=0;ix < this->pbasis;ix++) {
         h[ix] = std::real( h_prefactor[ix] * grad_rho[ix] );
     }

     for(int ig = 0;ig < this->pbasis;ig++) {
        double gradient2 = this->gx[ig]*this->gx[ig] + this->gy[ig]*this->gy[ig] + this->gz[ig]*this->gz[ig];
        if ( gradient2 > 0.0 ) h[ig] = h[ig] / sqrt( gradient2 );
     }

     pfft_execute_dft(plan_forward, (double (*)[2])h, (double (*)[2])h);

     for(int ix=0;ix < this->pbasis;ix++) {
         h[ix] = i * tpiba * this->plane_waves->g[ix].a[icar] * h[ix] * this->plane_waves->gmask[ix];
     }
//     if (gamma_only) h(nlm(:)) = CONJG(h(nl(:)))
     pfft_execute_dft(plan_back, (double (*)[2])h, (double (*)[2])h);
     double scale = 1.0 / (double)this->N;
     for(int ix=0;ix < this->pbasis;ix++) potential[ix] -= scale * std::real(h[ix]);

  }


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
      rFs = pow( 1 + 15.0*fa*(s*s) + fb*pow(s,4.0) + fc*pow(s,6.0) ,(1.0/15.0));
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
         rdFs_ds = ( 30.0*fa*s + 4.0*fb*pow(s,3.0) + 6.0*fc*pow(s,5.0) ) 
                / ( 15.0*pow(( 1.0 + 15.0*fa*(s*s) + fb*pow(s,4.0) + fc*pow(s, 6.0) ),(14.0/15.0)) );

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
    return pow(( 3.0 * PI*PI * rho ), (1.0/3.0));
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
       e = e + pow((q/q_cut),(double)ix)/(double)ix;
       dq0_dq = dq0_dq + pow((q/q_cut),(ix-1));
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
  if ( k >= (double)Vdw::Nrpoints*Vdw::dk ) {

      //std::cout << "GGGGGGG  " <<  k  << "  "  << (double)Vdw::Nrpoints*Vdw::dk << std::endl;
      throw RmgFatalException() << "k value requested is out of range in " << __FILE__ << " at line " << __LINE__ << "\n";

  }

  for(int ix = 0;ix < Vdw::Nqs*Vdw::Nqs;ix++) kernel_of_k[ix] = 0.0;

  // This integer division figures out which bin k is in since the kernel
  // is set on a uniform grid.
  int k_i = (int)(k / Vdw::dk);

  double dk2 = dk * dk;
  double A = (Vdw::dk * (k_i+1.0) - k) / Vdw::dk;
  double B = (k - Vdw::dk*k_i) / Vdw::dk;
  double C = (A*A*A - A) * dk / 6.0;
  double D = (B*B*B - B) * dk / 6.0;

  for(int q1_i = 0;q1_i < Vdw::Nqs;q1_i++) {
      for(int q2_i = 0;q2_i <= q1_i;q2_i++) {
      
          kernel_of_k[q1_i + q2_i*Vdw::Nqs] = A*kernel[k_i][q1_i][ q2_i] + B*kernel[k_i+1][q1_i][q2_i]
             +(C*d2phi_dk2[k_i][q1_i][q2_i] + D*d2phi_dk2[k_i+1][q1_i][q2_i]);

          kernel_of_k[q2_i + q1_i*Vdw::Nqs] = kernel_of_k[q1_i + q2_i*Vdw::Nqs];

      }
  }

}
