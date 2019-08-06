/*
 *
 * Copyright (c) 2019, Wenchang Lu
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

#ifndef RMG_LdaU_on_H
#define RMG_LdaU_on_H 1

#include "rmgtypedefs.h"
#include "LocalObject.h"
#include "BaseGrid.h"

class LdaU_on {

public:

    LdaU_on(LocalObject<double> &, BaseGrid &);
    ~LdaU_on(void);
    void calc_ns_occ(LocalObject<double> &LocalOrbital, double *mat_X, BaseGrid &);
    void app_vhubbard(LocalObject<double> &HL, BaseGrid &);
    void write_ldaU(void);

    LocalObject<double> *AtomicOrbital;
    int tot_orbitals_ldaU;
    int ldaU_m;
    double Ehub;
    double Ecorrect;

    double *ns_occ;    

    //  <AtomicOrbital | LocalOrbial>  global array
    double *Upsi_mat;

    double *Upsi_mat_local;

    std::vector<double> Hubbard_U;
    std::vector<double> Hubbard_J0;
    std::vector<double> Hubbard_alpha;
    std::vector<double> Hubbard_beta;

    double_2d_array Hubbard_J;

};

extern  LdaU_on *ldaU_on;
#endif
