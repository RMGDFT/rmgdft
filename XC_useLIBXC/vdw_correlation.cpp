/*

  This function calculates the non-local correlation contribution to the
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

#include "const.h"
#include "vdW.h"

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

const int Nqs = 20;
const double Vdw::q_mesh[20] = {
1.0e-5             , 0.0449420825586261, 0.0975593700991365, 0.159162633466142,
0.231286496836006, 0.315727667369529 , 0.414589693721418 , 0.530335368404141,
0.665848079422965, 0.824503639537924 , 1.010254382520950 , 1.227727621364570,
1.482340921174910, 1.780437058359530 , 2.129442028133640 , 2.538050036534580,
3.016440085356680, 3.576529545442460 , 4.232271035198720 , 5.0 };

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
  this->q_cut = q_mesh[Nqs];

  // Local storage
  this->total_rho = new double[pbasis];
  this->gx = new double[pbasis];
  this->gy = new double[pbasis];
  this->gz = new double[pbasis];
  this->q0 = new double[pbasis]();
  this->dq0_drho = new double[pbasis]();
  this->dq0_dgradrho = new double[pbasis]();
  this->thetas = new double[pbasis];


  
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

  // Initialize q0-related arrays.
  for(int ix = 0;ix < this->pbasis;ix++) {
      this->q0[ix] = this->q_cut;      
      this->dq0_drho[ix] = 0.0;
      this->dq0_dgradrho[ix] = 0.0;
  }

  for(int ix = 0;ix < this->pbasis;ix++) {

      double trho = this->total_rho[ix];
      

  }
  
}

