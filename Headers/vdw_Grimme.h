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

#ifndef RMG_vdw_Grimme_H
#define RMG_vdw_Grimme_H 1

namespace  vdw_Grimme
{
    double scale6 = 0.75;
    double damp = 20.0;
    double rcut = 200.0;

    std::unordered_map<std::string, const double> C6 ({
            {"H",      4.857},
            {"He",     2.775},
            {"Li",    55.853},
            {"Be",    55.853},
            {"B",    108.584},
            {"C",     60.710},
            {"N",     42.670},
            {"O",     24.284},
            {"F",     26.018},
            {"Ne",    21.855},
            {"Na",   198.087},
            {"Mg",   198.087},
            {"Al",   374.319},
            {"Si",   320.200},
            {"P",    271.980},
            {"S",    193.230},
            {"Cl",   175.885},
            {"Ar",   159.927},
            {"K",    374.666},
            {"Ca",   374.666},
            {"Sc",   374.666},
            {"Ti",   374.666},
            {"V",    374.666},
            {"Cr",   374.666},
            {"Mn",   374.666},
            {"Fe",   374.666},
            {"Co",   374.666},
            {"Ni",   374.666},
            {"Cu",   374.666},
            {"Zn",   374.666},
            {"Ga",   589.405},
            {"Ge",   593.221},
            {"As",   567.896},
            {"Se",   438.498},
            {"Br",   432.600},
            {"Kr",   416.642},
            {"Rb",   855.833},
            {"Sr",   855.833},
            {"Y",    855.833},
            {"Zr",   855.833},
            {"Nb",   855.833},
            {"Mo",   855.833},
            {"Tc",   855.833},
            {"Ru",   855.833},
            {"Rh",   855.833},
            {"Pd",   855.833},
            {"Ag",   855.833},
            {"Cd",   855.833},
            {"In",  1294.678},
            {"Sn",  1342.899},
            {"Sb",  1333.532},
            {"Te",  1101.101},
            {"I",   1092.775},
            {"Xe",  1040.391},
            {"Cs", 10937.246},
            {"Ba",  7874.678},
            {"La",  6114.381},
            {"Ce",  4880.348},
            {"Pr",  4880.348},
            {"Nd",  4880.348},
            {"Pm",  4880.348},
            {"Sm",  4880.348},
            {"Eu",  4880.348},
            {"Gd",  4880.348},
            {"Tb",  4880.348},
            {"Dy",  4880.348},
            {"Ho",  4880.348},
            {"Er",  4880.348},
            {"Tm",  4880.348},
            {"Yb",  4880.348},
            {"Lu",  4880.348},
            {"Hf",  3646.454},
            {"Ta",  2818.308},
            {"W",   2818.308},
            {"Re",  2818.308},
            {"Os",  2818.308},
            {"Ir",  2818.308},
            {"Pt",  2818.308},
            {"Au",  2818.308},
            {"Hg",  1990.022},
            {"Tl",  1986.206},
            {"Pb",  2191.161},
            {"Bi",  2204.274},
            {"Po",  1917.830},
            {"At",  1983.327},
            {"Rn",  1964.906},
    });

    std::unordered_map<std::string, const double> Ri ({
            {"H",   1.892},
            {"He",  1.912},
            {"Li",  1.559},
            {"Be",  2.661},
            {"B",   2.806},
            {"C",   2.744},
            {"N",   2.640},
            {"O",   2.536},
            {"F",   2.432},
            {"Ne",  2.349},
            {"Na",  2.162},
            {"Mg",  2.578},
            {"Al",  3.097},
            {"Si",  3.243},
            {"P",   3.222},
            {"S",   3.180},
            {"Cl",  3.097},
            {"Ar",  3.014},
            {"K",   2.806},
            {"Ca",  2.785},
            {"Sc",  2.952},
            {"Ti",  2.952},
            {"V",   2.952},
            {"Cr",  2.952},
            {"Mn",  2.952},
            {"Fe",  2.952},
            {"Co",  2.952},
            {"Ni",  2.952},
            {"Cu",  2.952},
            {"Zn",  2.952},
            {"Ga",  3.118},
            {"Ge",  3.264},
            {"As",  3.326},
            {"Se",  3.347},
            {"Br",  3.305},
            {"Kr",  3.264},
            {"Rb",  3.076},
            {"Sr",  3.035},
            {"Y",   3.097},
            {"Zr",  3.097},
            {"Nb",  3.097},
            {"Mo",  3.097},
            {"Tc",  3.097},
            {"Ru",  3.097},
            {"Rh",  3.097},
            {"Pd",  3.097},
            {"Ag",  3.097},
            {"Cd",  3.097},
            {"In",  3.160},
            {"Sn",  3.409},
            {"Sb",  3.555},
            {"Te",  3.575},
            {"I",   3.575},
            {"Xe",  3.555},
            {"Cs",  3.405},
            {"Ba",  3.330},
            {"La",  3.251},
            {"Ce",  3.313},
            {"Pr",  3.313},
            {"Nd",  3.313},
            {"Pm",  3.313},
            {"Sm",  3.313},
            {"Eu",  3.313},
            {"Gd",  3.313},
            {"Tb",  3.313},
            {"Dy",  3.313},
            {"Ho",  3.313},
            {"Er",  3.313},
            {"Tm",  3.313},
            {"Yb",  3.313},
            {"Lu",  3.313},
            {"Hf",  3.378},
            {"Ta",  3.349},
            {"W",   3.349},
            {"Re",  3.349},
            {"Os",  3.349},
            {"Ir",  3.349},
            {"Pt",  3.349},
            {"Au",  3.349},
            {"Hg",  3.322},
            {"Tl",  3.752},
            {"Pb",  3.673},
            {"Bi",  3.586},
            {"Po",  3.789},
            {"At",  3.762},
            {"Rn",  3.636},
    });

}

#endif
