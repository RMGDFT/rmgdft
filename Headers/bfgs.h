/*

			bfgs.h

  Interface to QE fortran 90 bfgs module.
*/

#include "rmgtypedefs.h"
#include "rmg_mangling.h"

#define init_bfgs       RMG_FC_MODULE(bfgs_module,init_bfgs,mod_BFGS_MODULE,INIT_BFGS)
#define	bfgs    	RMG_FC_MODULE(bfgs_module,bfgs,mod_BFGS_MODULE,BFGS)
#define	reset_bfgs    	RMG_FC_MODULE(bfgs_module,reset_bfgs,mod_BFGS_MODULE,RESET_BFGS)
#define	terminate_bfgs 	RMG_FC_MODULE(bfgs_module,terminate_bfgs,mod_BFGS_MODULE,TERMINATE_BFGS)

#if __cplusplus
extern "C" {
#endif

void init_bfgs( int *stdout_, int *bfgs_ndim_, double *trust_radius_max_,
               double *trust_radius_min_, double *trust_radius_ini_, 
               double *w_1_, double *w_2_);

// h, fcell and iforce_h are 3x3 double matrices.
// pos_in and grad_in are linear arrays of doubles sized a bit larger than
// 3*nat
//
//void bfgs( filebfgs, double *pos_in, h, double *nelec, double *energy,
void bfgs( double *pos_in, double *h, double *nelec, double *energy,
           double *grad_in, double *fcell, double *iforceh, double *felec,
           double *energy_thr, double *grad_thr, double *cell_thr, double *fcp_thr,
           double *energy_error, double *grad_error, double *cell_error, double *fcp_error,
           bool *lmovecell, bool *lfcp, double *fcp_cap, double *fcp_hess, bool *step_accepted,
           bool *stop_bfgs, bool *failed, int *istep );

void reset_bfgs( int *, bool *lfcp, double *fcp_hess );

void terminate_bfgs( double *energy, double *energy_thr, double *grad_thr, 
                     double *cell_thr, double *fcp_thr,
                     bool *lmovecell, bool *lfcp );

#if __cplusplus
}
#endif
