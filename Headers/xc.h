#ifndef RMG_XC_H
#define RMG_XC_H 1


void corlyp (rmg_double_t * dp, rmg_double_t * dm, rmg_double_t * dp1, rmg_double_t * dm1, rmg_double_t * dp2, rmg_double_t * dm2, rmg_double_t * ec, rmg_double_t * vcp0, rmg_double_t * vcm0, int *ndm);
void xclda_pz81 (rmg_double_t * rho, rmg_double_t * vxc_f, int n);
void exclda_pz81 (rmg_double_t * rho, rmg_double_t * exc, int n);


#endif
