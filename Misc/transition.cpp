// Some transitional routines


#include "rmgtypes.h"
#include "mpi.h"
#include "transition.h"
#include "params.h"
#include "rmgtypes.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Subdiag.h"

BaseGrid *Rmg_G;
TradeImages *Rmg_T;
Lattice Rmg_L;

// Pointer to Kpoint class array
/* Hartree potential */
extern double *vh;
/* Nuclear local potential */
extern double *vnuc;
/* Exchange-correlation potential */
extern double *vxc;
// Pointer to Kpoint class array
extern void **Kptr;
/*  Electronic charge density of pposite spin density*/
extern double *rho_oppo;
/* Core Charge density */
extern double *rhocore;
/* Compensating charge density */
extern double *rhoc;



extern "C"
MPI_Comm transition_get_grid_comm(void)
{
    return pct.grid_comm;
}
extern "C" int transition_get_gridpe(void)
{
    return pct.gridpe;
}

extern "C" void subdiag_gamma (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc)
{
     Kpoint<double>* kptr = (Kpoint<double> *)Kptr[0]; 
     Subdiag<double> (kptr, vh, vnuc, vxc, ct.subdiag_driver);
}

extern "C" bool scf (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc,
          rmg_double_t * rho, rmg_double_t * rho_oppo, rmg_double_t * rhocore, rmg_double_t * rhoc )
{

     if(ct.is_gamma) {
         Kpoint<double> **kptr = (Kpoint<double> **)Kptr; 
         return Scf<double> (vxc, vh, ct.vh_ext,
              vnuc, rho, rho_oppo, rhocore, rhoc, ct.spin_flag,
              ct.hartree_min_sweeps, ct.hartree_max_sweeps , ct.boundaryflag, kptr);
     }
     else {
         Kpoint<std::complex<double> > **kptr = (Kpoint<std::complex<double> > **)Kptr; 
         return Scf<std::complex<double> >(vxc, vh, ct.vh_ext,
              vnuc, rho, rho_oppo, rhocore, rhoc, ct.spin_flag,
              ct.hartree_min_sweeps, ct.hartree_max_sweeps , ct.boundaryflag, kptr);
     }
}

extern "C" void get_new_rho (STATE * states, rmg_double_t * rho)
{

     if(ct.is_gamma) {
         Kpoint<double> **kptr = (Kpoint<double> **)Kptr; 
         GetNewRho<double> (kptr, rho);
     }
     else {
         Kpoint<std::complex<double> > **kptr = (Kpoint<std::complex<double> >**)Kptr;
         GetNewRho<std::complex<double> >(kptr, rho);
     }

}

extern "C" void mg_eig_state(STATE *sp, int tid, double *vtot_psi)
{
    Kpoint<double> *kptr = (Kpoint<double> *)Kptr[0];
    MgEigState<double,float> (kptr, sp, tid, vtot_psi);
}

extern "C" void mg_eig_state_driver (STATE * sp, int tid, double * vtot_psi)
{

        if(ct.is_gamma) {
            Kpoint<double> *kptr = (Kpoint<double> *)Kptr[0];
            MgEigState<double,float> (kptr, sp, tid, vtot_psi);
        }
        else {
            Kpoint<std::complex<double>> *kptr = (Kpoint<std::complex<double>> *)Kptr[0];
            MgEigState<std::complex<double>,std::complex<float> > (kptr, sp, tid, vtot_psi);
        }


}

