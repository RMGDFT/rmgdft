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
#include "const.h"
#include "vdW.h"
#include "xc.h"

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

int Vdw::Nqs = 20;
double Vdw::q_mesh[20] = {
1.0e-5             , 0.0449420825586261, 0.0975593700991365, 0.159162633466142,
0.231286496836006, 0.315727667369529 , 0.414589693721418 , 0.530335368404141,
0.665848079422965, 0.824503639537924 , 1.010254382520950 , 1.227727621364570,
1.482340921174910, 1.780437058359530 , 2.129442028133640 , 2.538050036534580,
3.016440085356680, 3.576529545442460 , 4.232271035198720 , 5.0 };
const double Vdw::epsr = 1.0e-12;

/*

  rho_valence - valence charge
  rho_core    - core correction if applicable
  

*/

Vdw::Vdw (BaseGrid &G, Lattice &L, TradeImages &T, int type, double *rho_valence, double *rho_core, double &etxc, double &vtxc, double *v)
{

  // Grid parameters
  this->pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
  this->hxgrid = G.get_hxgrid(G.default_FG_RATIO);
  this->hygrid = G.get_hygrid(G.default_FG_RATIO);
  this->hzgrid = G.get_hzgrid(G.default_FG_RATIO);
  this->dimx = G.get_PX0_GRID(G.default_FG_RATIO);
  this->dimy = G.get_PY0_GRID(G.default_FG_RATIO);
  this->dimz = G.get_PZ0_GRID(G.default_FG_RATIO);

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
  this->thetas = new std::complex<double> [this->pbasis*Vdw::Nqs];


  
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

  get_q0_on_grid ();

}


// Destructor just frees memory
Vdw::~Vdw(void)
{

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
          printf("THETAS for iq=%d  ix=%d is %14.12f\n",iq,ix,thetas[ix + iq*this->pbasis]);
      }
  }

//  do idx = 1, Nqs
//     CALL fwfft ('Dense', thetas(:,idx), dfftp)
//  end do

  
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
      rFs = pow( 1 + 15.0*fa*(s*s) + fb*(s*s*s*s) + fc*(s*s*s*s*s*s) ,(1.0/15.0));
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
         rdFs_ds = ( 30.0*fa*s + 4.0*fb*(s*s*s) + 6.0*fc*(s*s*s*s*s) ) 
                / ( 15.0*pow(( 1.0 + 15.0*fa*(s*s) + fb*(s*s*s*s) + fc*(s*s*s*s*s*s) ),(14.0/15.0)) );

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

