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

// Some transitional routines


#include "mpi.h"
#include "transition.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Subdiag.h"

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

#if 0
extern "C" void subdiag_gamma (STATE * states, double * vh, double * vnuc, double * vxc)
{
     Kpoint<double>* kptr = (Kpoint<double> *)Kptr[0]; 
     Subdiag<double> (kptr, vh, vnuc, vxc, ct.subdiag_driver);
}

extern "C" bool scf (STATE * states, double * vxc, double * vh, double * vnuc,
          double * rho, double * rho_oppo, double * rhocore, double * rhoc )
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

extern "C" void get_new_rho (STATE * states, double * rho)
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
    MgEigState<double,float> (kptr, sp, vtot_psi);
}

extern "C" void mg_eig_state_driver (STATE * sp, int tid, double * vtot_psi)
{

        if(ct.is_gamma) {
            Kpoint<double> *kptr = (Kpoint<double> *)Kptr[0];
            MgEigState<double,float> (kptr, sp, vtot_psi);
        }
        else {
            Kpoint<std::complex<double>> *kptr = (Kpoint<std::complex<double>> *)Kptr[0];
            MgEigState<std::complex<double>,std::complex<float> > (kptr, sp, vtot_psi);
        }


}
#endif
