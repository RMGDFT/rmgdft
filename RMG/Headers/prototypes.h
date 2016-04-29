/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include <stdbool.h>

#if __cplusplus
extern "C" {
#endif

void CompatRmgTimerPrint(const char *outfile, int steps);
void betaxpsi1_calculate_one(STATE *st, int ion, int nion, double *sintR, double *sintI, int kpt, double *weiptr_base);
void get_te (double *rho, double *rho_oppo, double *rhocore, double *rhoc, double *vh, double *vxc,
             STATE *states, int ii_flag);
void get_vxc (double *rho, double *rho_oppo, double *rhocore, double *vxc);
void nlforce_par_gamma (double * par_gamma, int ion, int nh, double *force);
void nlforce_par_omega (double * par_omega, int ion, int nh, double *force);
void nlforce_par_Q (double *veff, double *, double *, double *gamma, int ion, ION *iptr, int nh,
                     double *forces);
void subdiag_app_B_one (STATE *sp, double * b_psi);
void subdiag_app_A_one (STATE *sp, double * a_psi, double * s_psi, double * vtot_eig);
void subdiag_app_AB_one (STATE *sp, double * a_psi, double * b_psi, double * vtot_eig_s);
void subdiag_app_AB (STATE * states, double * a_psi, double * b_psi, double * vtot_eig);
#if __cplusplus
}
#endif

