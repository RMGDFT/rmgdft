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

#ifndef RMG_State_H
#define RMG_State_H 1

#include <mpi.h>
template <typename T> class Kpoint;
template <typename StateType> class State {

public:
    State(void);

    void set_storage(StateType *storage);
    void normalize(StateType *tpsi, int istate);
    bool is_occupied(void);

    // kpoint this state is attached to
    Kpoint<StateType> *Kptr;

    /** Index showing which k-point this orbital is associated with */
    int kidx;

    //void normalize(void);

    // Storage area for the orbital
    StateType *psi;
 
    /** Nuclear potential */
    double *vnuc;

    /** Hartree potential */
    double *vh;

    /** Exchange correlation potential */
    double *vxc;

    /** Total potential */
    double *vtot;

    /** Total basis size on each processor (dimx*dimy*dimz) */
    int pbasis;

    /** Wavefunction residual error computed by multigrid solver */
    double res;

    // Last two eigenvalues
    double eig[2];
    double oldeig[2];

    // First computed eigenvalue from MgEigState per scf cycle
    double feig[2];

    // Index of the orbital 
    int istate;

    // Occupation of the orbital
    double occupation[2];


};

#endif
