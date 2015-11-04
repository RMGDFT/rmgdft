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
#include "Lattice.h"
#include "TradeImages.h"
#include "FiniteDiff.h"
#include "vdW.h"

/*

  rho_valence - valence charge
  rho_core    - core correction if applicable
  

*/

void xc_vdW_DF (BaseGrid &G, Lattice &L, TradeImages &T, int type, double *rho_valence, double *rho_core, double &etxc, double &vtxc, double *v)
{

  int pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
  double hxgrid = G.get_hxgrid(G.default_FG_RATIO);
  double hygrid = G.get_hygrid(G.default_FG_RATIO);
  double hzgrid = G.get_hzgrid(G.default_FG_RATIO);
  int dimx = G.get_PX0_GRID(G.default_FG_RATIO);
  int dimy = G.get_PY0_GRID(G.default_FG_RATIO);
  int dimz = G.get_PZ0_GRID(G.default_FG_RATIO);



  // Local storage
  double *total_rho = new double[pbasis];
  double *gx = new double[pbasis];
  double *gy = new double[pbasis];
  double *gz = new double[pbasis];


  
  // Get total charge and compute it's gradient
  for(int i = 0;i < pbasis;i++) total_rho[i] = rho_valence[i] + rho_core[i];

  CPP_app_grad_driver (&L, &T, total_rho, gx, gy, gz, dimx, dimy, dimz, hxgrid, hygrid, hzgrid, APP_CI_SIXTH);
  



  delete [] gz;
  delete [] gy;
  delete [] gx;
  delete [] total_rho;

}
