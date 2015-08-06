#ifndef RMG_XC_H
#define RMG_XC_H 1


void corlyp (double * dp, double * dm, double * dp1, double * dm1, double * dp2, double * dm2, double * ec, double * vcp0, double * vcm0, int *ndm);
void xclda_pz81 (double * rho, double * vxc_f, int n);
void exclda_pz81 (double * rho, double * exc, int n);


#endif
