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

#ifndef RMG_Prolong_H
#define RMG_Prolong_H 1

#include <vector>
#include "TradeImages.h"
#include "rmg_error.h"

#define MAX_PROLONG_ORDER 10

class Prolong {

public:
    Prolong(int ratio, int order, TradeImages &TI);
    ~Prolong(void);
    template<typename T>
    void prolong (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz);
    template<typename T>
    void prolong_ortho_10 (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz);


private:
    void cgen_prolong (double *coef, double fraction);

    int ratio;
    int order;
    TradeImages &TR;
    double a[MAX_PROLONG_ORDER][MAX_PROLONG_ORDER];
};

#endif
