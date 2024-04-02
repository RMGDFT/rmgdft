#ifndef RMG_qk_compat_H
#define RMG_qk_compat_H 1


#include "rmg_mangling.h"

#define kpoint_grid  RMG_FC_GLOBAL(kpoint_grid, KPOINT_GRID)
#define irreducible_bz  RMG_FC_GLOBAL(irreducible_bz, IRREDUCIBLE_BZ)


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
