/*
 *
 * Copyright (c) 2019, Emil Briggs
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

#ifndef RMG_LdaU_H
#define RMG_LdaU_H 1

#include "rmgtypedefs.h"
#include "Kpoint.h"

template <typename KpointType> class LdaU {

public:

    LdaU(Kpoint<KpointType> &kp);
    ~LdaU(void);
    void calc_ns_occ(KpointType *sint, KpointType *sint_in);
    void app_vhubbard(KpointType *v_hub_x_psi, KpointType *sint);
    void write_ldaU(void);

    int ldaU_m;
    double Ehub;
    double Ecorrect;

    Kpoint<KpointType> &K;
    doubleC_4d_array ns_occ;    

    std::vector<double> Hubbard_U;
    std::vector<double> Hubbard_J0;
    std::vector<double> Hubbard_alpha;
    std::vector<double> Hubbard_beta;

    double_2d_array Hubbard_J;

};

#endif
