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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "blas.h"
#include <complex>
#include "transition.h"
#include "../Headers/prototypes.h"


template void Betaxpsi<double> (BaseGrid *, TradeImages *, Lattice *, Kpoint<double> *);
template void Betaxpsi<std::complex<double> > (BaseGrid *, TradeImages *, Lattice *, Kpoint<std::complex<double>> *);

// Transitional functions used to pack orbital array
void pack_to_complex(double *psi, int nstates, int pbasis)
{
    int offset=0;

    double *tarr = new double[2*pbasis];
    for(int st = 0;st < nstates;st++) {

        for(int idx = 0;idx < 2*pbasis;idx++) {
            tarr[idx] = psi[idx + offset];
        }

        for(int idx = 0;idx < pbasis;idx++) {
            psi[2*idx + offset] = tarr[idx];
            psi[2*idx + offset + 1] = tarr[idx + pbasis];
        }

        offset += 2*pbasis;
    }

    delete [] tarr;
}

void pack_to_standard(double *psi, int nstates, int pbasis)
{
    int offset=0;

    double *tarr = new double[2*pbasis];
    for(int st = 0;st < nstates;st++) {

        for(int idx = 0;idx < 2*pbasis;idx++) {
            tarr[idx] = psi[idx + offset];
        }

        for(int idx = 0;idx < pbasis;idx++) {
            psi[idx + offset] = tarr[2*idx];
            psi[idx + pbasis + offset] = tarr[2*idx + 1];
        }

        offset += 2*pbasis;
    }

    delete [] tarr;
}


// This is currently a transitional implementation. It transforms the wavefunction arrays
// to the old storage format, applies the operators and then transforms the output array
// to the new storage format.
template <typename OrbitalType>
void Betaxpsi (BaseGrid *G, TradeImages *T, Lattice *L, Kpoint<OrbitalType> *Kptr)
{

    if(typeid(OrbitalType) == typeid(std::complex<double>)) {
        pack_to_standard((double *)Kptr->orbital_storage, Kptr->nstates, Kptr->pbasis);
        betaxpsi1 (Kptr->kstates, Kptr->kidx);
        pack_to_complex((double *)Kptr->orbital_storage, Kptr->nstates, Kptr->pbasis);
    }
    else {
        betaxpsi1 (Kptr->kstates, Kptr->kidx);
    }

}

