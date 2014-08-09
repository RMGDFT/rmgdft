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

#include "transition.h"

void MixBetaxpsi (int mix, int kpt)
{

    int size = pct.num_nonloc_ions * ct.num_states * ct.max_nl;
    int offset = kpt * size;

    double *newsintR = &pct.newsintR_local[offset];
    double *oldsintR = &pct.oldsintR_local[offset];
    double *newsintI = &pct.newsintI_local[offset];
    double *oldsintI = &pct.oldsintI_local[offset];
    double scale = 1.0 - ct.prjmix;

    if (mix)
    {
        for(int idx = 0;idx < size;idx++) {
            oldsintR[idx] = scale * oldsintR[idx];
            oldsintR[idx] = oldsintR[idx] + ct.prjmix * newsintR[idx];
        }

        if(!ct.is_gamma) {
            for(int idx = 0;idx < size;idx++) {
                oldsintI[idx] = scale * oldsintI[idx];
                oldsintI[idx] = oldsintI[idx] + ct.prjmix * newsintI[idx];
            }
        }
    }


    else
    {
        for(int idx = 0;idx < size;idx++)
            oldsintR[idx] = newsintR[idx];

        if(!ct.is_gamma) {
            for(int idx = 0;idx < size;idx++)
                oldsintI[idx] = newsintI[idx];
        }
    }

}

