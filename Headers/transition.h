#ifndef RMG_transition_h
#define RMG_transition_h

#if __cplusplus
#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"


extern BaseGrid *Rmg_G;
extern TradeImages *Rmg_T;
extern Lattice Rmg_L;

extern "C"
{
double my_crtc (void);
MPI_Comm transition_get_grid_comm(void);
void thread_barrier_wait(void);
int transition_get_gridpe(void);
void get_vxc (double * rho, double * rho_oppo, double * rhocore, double * vxc);
}
extern "C" void get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target, int boundaryflag);
extern "C" void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio);
extern "C" void mix_rho (double * new_rho, double * rho, double *rhocore, int length, int length_x, int length_y, int length_z);
extern "C" void  get_rho_oppo (double * rho, double * rho_oppo);
extern "C" void get_ddd (double *veff);
extern "C" void mix_betaxpsi (int mix);
extern "C" void rmg_lbfgs (void);
extern "C" void write_restart (char *name, double * vh, double * rho, double * rho_oppo, double * vxc, STATE * states);

template <typename OrbitalType> void Init (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> void Relax (int steps, double * vxc, double * vh, double * vnuc,
              double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> bool Quench (double * vxc, double * vh, double * vnuc, double * rho,
             double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> bool Scf (double * vxc, double * vh, double *vh_ext,
          double * vnuc, double * rho, double * rho_oppo, double * rhocore, double * rhoc, int spin_flag,
          int hartree_min_sweeps, int hartree_max_sweeps , int boundaryflag, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> void AppNls(Kpoint<OrbitalType> *kpoint, double *sintR, double *sintI);
template <typename OrbitalType, typename CalcType> void MgEigState (BaseGrid *G, TradeImages *T, Lattice *L, STATE * sp, int tid, double * vtot_psi);




#endif
#endif

#define rmg_printf( message... ) \
         fprintf( ct.logfile, message )

