/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#ifndef RMG_Solvers_H
#define RMG_Solvers_H 1

#if __cplusplus
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "Kpoint.h"

template <typename OrbitalType> void MgridSubspace (Kpoint<OrbitalType> *kptr, double *vtot);
template <typename OrbitalType> void Davidson(Kpoint<OrbitalType> *kptr, double *vtot, int &notconv);


template <typename KpointType>
double ApplyHamiltonian (Kpoint<KpointType> *kptr, KpointType *psi, KpointType *h_psi, double *vtot, KpointType *nv);
template <typename KpointType>
double ApplyHamiltonianBlock (Kpoint<KpointType> *kptr, int first_state, int num_states, KpointType *h_psi, double *vtot);
template <typename RmgType>
void CPP_app_smooth_test (RmgType * f, RmgType *b, int dimx, int dimy, int dimz);
template <typename OrbitalType>
void DavPreconditioner (Kpoint<OrbitalType> *kptr, OrbitalType *res, 
                        double fd_diag, double *eigs, double *vtot, int notconv, double avg_potential);
template <typename OrbitalType>
void DavPreconditionerOne (Kpoint<OrbitalType> *kptr, OrbitalType *res, 
                        double fd_diag, double eig, double *vtot, double avg_potential);




#endif
#endif


