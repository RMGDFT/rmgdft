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
    void calc_ns_occ(KpointType *sint, int first_state, int nstates);
    void calc_force(KpointType *sint, double *force_ldau);
    void app_vhubbard(KpointType *v_hub_x_psi, KpointType *sint, int first_state, int nstates);
    void write_ldaU(void);
    void init_ns_occ(void);
    void Hubbard_matrix();

    void calc_energy();
    int ldaU_m;
    double Ehub;
    double Ecorrect;

    Kpoint<KpointType> &K;
    boost::multi_array<std::complex<double>, 4> ns_occ;
    boost::multi_array<std::complex<double>, 4> vhub_ns;

    std::vector<double> Hubbard_U;
    std::vector<double> Hubbard_J0;
    std::vector<double> Hubbard_alpha;
    std::vector<double> Hubbard_beta;

    double_2d_array Hubbard_J;

    std::unordered_map<std::string, const double> hubbard_occ={
            {"H",   1.0},
            {"He",  0.0},
            {"Li",  0.0},
            {"Be",  0.0},
            {"B",   0.0},
            {"C",   2.0},
            {"N",   3.0},
            {"O",   4.0},
            {"F",   0.0},
            {"Ne",  0.0},
            {"Na",  0.0},
            {"Mg",  0.0},
            {"Al",  0.0},
            {"Si",  0.0},
            {"P",   0.0},
            {"S",   0.0},
            {"Cl",  0.0},
            {"Ar",  0.0},
            {"K",   0.0},
            {"Ca",  0.0},
            {"Sc",  0.0},
            {"Ti",  2.0},
            {"V",   3.0},
            {"Cr",  4.0},
            {"Mn",  5.0},
            {"Fe",  6.0},
            {"Co",  7.0},
            {"Ni",  8.0},
            {"Cu", 10.0},
            {"Zn", 10.0},
            {"Ga", 10.0},
            {"Ge",  0.0},
            {"As",  0.0},
            {"Se",  0.0},
            {"Br",  0.0},
            {"Kr",  0.0},
            {"Rb",  0.0},
            {"Sr",  0.0},
            {"Y",   0.0},
            {"Zr",  2.0},
            {"Nb",  3.0},
            {"Mo",  4.0},
            {"Tc",  5.0},
            {"Ru",  6.0},
            {"Rh",  7.0},
            {"Pd",  8.0},
            {"Ag", 10.0},
            {"Cd", 10.0},
            {"In", 10.0},
            {"Sn",  0.0},
            {"Sb",  0.0},
            {"Te",  0.0},
            {"I",   0.0},
            {"Xe",  0.0},
            {"Cs",  0.0},
            {"Ba",  0.0},
            {"La",  0.0},
            {"Ce",  2.0},
            {"Pr",  3.0},
            {"Nd",  4.0},
            {"Pm",  5.0},
            {"Sm",  6.0},
            {"Eu",  6.0},
            {"Gd",  7.0},
            {"Tb",  8.0},
            {"Dy",  9.0},
            {"Ho", 10.0},
            {"Er", 11.0},
            {"Tm", 12.0},
            {"Yb", 13.0},
            {"Lu", 14.0},
            {"Hf",  2.0},
            {"Ta",  3.0},
            {"W",   4.0},
            {"Re",  5.0},
            {"Os",  6.0},
            {"Ir",  7.0},
            {"Pt",  8.0},
            {"Au", 10.0},
            {"Hg", 10.0},
            {"Tl",  0.0},
            {"Pb",  0.0},
            {"Bi",  0.0},
            {"Po",  0.0},
            {"At",  0.0},
            {"Rn",  0.0},
            {"Th",  2.0},
            {"Pa",  3.0},
            {"U",  4.0},
            {"Np",  5.0},
            {"Pu",  6.0},
            {"Am",  6.0},
            {"Cm",  7.0},
            {"Bk",  8.0},
            {"Cf",  9.0},
            {"Es", 10.0},
            {"Fm", 11.0},
            {"Md", 12.0},
            {"No", 13.0},
            {"Lr", 14.0},
    };

};

#endif
