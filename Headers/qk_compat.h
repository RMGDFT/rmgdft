#ifndef RMG_qk_compat_H
#define RMG_qk_compat_H 1


#include "rmg_mangling.h"


#define inverse_s	RMG_FC_MODULE(symm_base,inverse_s,mod_SYMM_BASE,INVERSE_S)
#define set_sym_bl	RMG_FC_MODULE(symm_base,set_sym_bl,mod_SYMM_BASE,SET_SYM_BL)
#define find_sym	RMG_FC_MODULE(symm_base,find_sym,mod_SYMM_BASE,FIND_SYM)
#define sgam_at	        RMG_FC_MODULE(symm_base,sgam_at,mod_SYMM_BASE,SGAM_AT)
#define sgam_at_mag	RMG_FC_MODULE(symm_base,sgam_at_mag,mod_SYMM_BASE,SGAM_AT_MAG)
#define set_sym	        RMG_FC_MODULE(symm_base,set_sym,mod_SYMM_BASE,SET_SYM)
#define checkallsym	RMG_FC_MODULE(symm_base,checkallsym,mod_SYMM_BASE,CHECKALLSYM)
#define s_axis_to_cart	RMG_FC_MODULE(symm_base,s_axis_to_cart,mod_SYMM_BASE,S_AXIS_TO_CART)
#define find_sym_ifc	RMG_FC_MODULE(symm_base,find_sym_ifc,mod_SYMM_BASE,FIND_SYM_IFC)
#define sgam_at_ifc	RMG_FC_MODULE(symm_base,sgam_at_ifc,mod_SYMM_BASE,SGAM_AT_IFC)
#define remove_sym	RMG_FC_MODULE(symm_base,remove_sym,mod_SYMM_BASE,REMOVE_SYM)
#define kpoint_grid     RMG_FC_GLOBAL(kpoint_grid, KPOINT_GRID)
#define irreducible_bz  RMG_FC_GLOBAL(irreducible_bz, IRREDUCIBLE_BZ)

extern "C" void set_sym_bl(double *at, double *bg, int *nrot_out);

extern "C" void find_sym( int *nat, double *tau, int *ityp, int *magnetic_sym, double *m_loc, int *no_z_inv );

extern "C" void kpoint_grid( int *nrot, int *time_reversal, int *skip_equivalence,
                  int *s, int *t_rev, double *, int *npk,
                  int *k1, int *k2, int *k3, int *nk1,int *nk2,int *nk3,
                  int *nks, double *xk, double *wk );

extern "C" void irreducible_bz( int *nrot,
                           int *s,
                           int *nsym,
                           int *minus_q,
                           int *magnetic_sym,
                           double *at,         // Use Rmg_L.at
                           double *bg,         // Use Rmg_L.bg
                           int *npk,           // Max number of k-points
                           int *nks,           // Number of k-points
                           double *xk,         // K-points stored as triplets
                           double *wk,         // Their weights
                           int *t_rev );       // Use Rmg_Symm.time_rev.data()

#endif
