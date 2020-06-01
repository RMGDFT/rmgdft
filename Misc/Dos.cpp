/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
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

#include <cstdint>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "GlobalSums.h"
#include "Symmetry.h"
#include "rmg_error.h"
#include "Dos.h"
#include "transition.h"


Dos::Dos (int kmesh[3], int kshift[3], Lattice &L_in, double gaus_broad_in) : L(L_in), gaus_broad(gaus_broad_in)
{

    tetraflag = false;
    if(ct.dos_flag == 0)  tetraflag = true;
    if(kmesh[0] * kmesh[1] * kmesh[2] < 1) tetraflag = false;  

    if(tetraflag) 
    {
        set_Pij();
        tetra_init(kmesh, kshift);
    }

}
void Dos::tot_dos(int nk, int nband, std::vector<double> eigs, double Ef)
{
    double Emax = *std::max_element(eigs.begin(), eigs.end());
    double Emin = *std::min_element(eigs.begin(), eigs.end());
    Emin = std::floor(Emin - Ef);
    Emax = std::ceil(Emax - Ef);
    std::cout << "gaus " << gaus_broad << std::endl;
    double delta_e = 0.01;
    int num_epoint = (int) ((Emax - Emin)/delta_e) + 1;
    std::vector<double> dos_t, e_list;
    dos_t.resize(num_epoint);
    e_list.resize(num_epoint);
    std::cout << " num " << num_epoint<< std::endl;
    for(int ie = 0; ie < num_epoint; ie++)
    {
        e_list[ie] = Emin + ie * delta_e;
        double energy = Emin + ie * delta_e + Ef;
        if(tetraflag)
            dos_t[ie] = tot_dos_tetra(nk, nband, eigs, energy);
        else
            dos_t[ie] = tot_dos_gauss(nk, nband, eigs, energy);
    }

    if(pct.gridpe == 0)
    {
        std::string filename = "dos_tot_spin" + std::to_string(pct.spinpe);
        FILE *dos_f = fopen (filename.c_str(), "w");
        for(int ie = 0; ie < num_epoint; ie++)
        {
            fprintf(dos_f, " %e  %e \n", e_list[ie], dos_t[ie]);
        }
        fclose(dos_f);
    }

}
double Dos::tot_dos_gauss(int nk, int nband, std::vector<double> eigs, double energy)
{
    double dos_one = 0.0;
    if(gaus_broad < 1.0e-5) gaus_broad = 0.1;
    if(eigs.size() < (size_t)nk * nband)
    {
        throw RmgFatalException() << "eigs vector is not large enough" << " in " << __FILE__ << " at line " << __LINE__ << ".\n";
    }
    for(int ik = 0; ik < nk; ik++)
    {
        for(int ib = 0; ib < nband; ib++)
        {
            double tem = energy - eigs[ik * nband + ib];
            dos_one += std::exp(-tem * tem/gaus_broad/gaus_broad) / gaus_broad /std::sqrt(PI) * ct.kp[ik].kweight;
        }
    }
    return dos_one;
}
double Dos::tot_dos_tetra(int nk, int nband, std::vector<double> eigs, double energy)
{
    double dos_one = 0.0;
    std::vector<double> e_tetra;
    if(eigs.size() < (size_t)nk * nband)
    {
        throw RmgFatalException() << "eigs vector is not large enough" << " in " << __FILE__ << " at line " << __LINE__ << ".\n";
    }

    double e1, e2, e3, e4;
    for(int it = 0; it < num_tetra; it++)
    {
        for(int ib = 0; ib < nband; ib++)
        {
            e_tetra.assign(4, 0.0);

            for(int ip = 0; ip < nk_per_tetra; ip++)
            {
                int ik = tetra[ip][it];
                for(int id = 0; id < 4; id++)
                    e_tetra[id] += Pij[id][ip] * eigs[ik * nband + ib];
            }
            std::sort(e_tetra.begin(), e_tetra.end());
            e1 = e_tetra[0];
            e2 = e_tetra[1];
            e3 = e_tetra[2];
            e4 = e_tetra[3];
            if(energy >= e4 || energy <= e1 ) continue;
            if(energy >= e3)
            {
                dos_one += 3.0 * (e4-energy) *  (e4-energy) /(e4-e1)/(e4-e2)/(e4-e3)/ num_tetra;
            }
            else if(energy >= e2)
            {
                double tem = 3.0 *(e2-e1) + 6.0 *(energy - e2) 
                    -3.0 * (e3-e1+e4-e2)/(e3-e2)/(e4-e2) *(energy-e2) * (energy-e2);
                dos_one += tem/(e3-e1)/(e4-e1)/num_tetra;
            }
            else 
            {
                dos_one += 3.0 * (energy-e1) *  (energy-e1) /(e4-e1)/(e3-e1)/(e2-e1)/ num_tetra;
            }

        }
    }
    if(ct.nspin == 1) dos_one *=2.0;
    return dos_one;

}
void Dos::tetra_init(int kmesh[3], int kshift[3])
{
    std::vector<int> kequiv;
    kequiv.resize(kmesh[0] * kmesh[1] * kmesh[2]);
    num_tetra =kmesh[0] * kmesh[1] * kmesh[2] * 6;
    nk_per_tetra = 20;
    tetra.resize(boost::extents[nk_per_tetra][num_tetra]);

    // find map between all kpoints and irreducible kpoints

    double sym_qvec[3], dk[3], xk[3];
    for(int i = 0; i < kmesh[0]; i++)
    {
        for(int j = 0; j < kmesh[1]; j++)
        {
            for(int k = 0; k < kmesh[2]; k++)
            {
                xk[0] = (i+ 0.5 * kshift[0])/kmesh[0];
                xk[1] = (j+ 0.5 * kshift[1])/kmesh[1];
                xk[2] = (k+ 0.5 * kshift[2])/kmesh[2];

                bool find_eq (false);
                for(int ik = 0; ik < ct.num_kpts; ik++)
                {
                    for(int isym = 0; isym < Rmg_Symm->nsym; isym++)
                    {
                        sym_qvec[0] = Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                            Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                            Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 2 ] * ct.kp[ik].kpt[2];
                        sym_qvec[1] = Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                            Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                            Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 2 ] * ct.kp[ik].kpt[2];
                        sym_qvec[2] = Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                            Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                            Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 2 ] * ct.kp[ik].kpt[2];

                        if(Rmg_Symm->time_rev[isym])
                        {
                            sym_qvec[0] *= -1.0;
                            sym_qvec[1] *= -1.0;
                            sym_qvec[2] *= -1.0;
                        }

                        dk[0] = sym_qvec[0] - xk[0];
                        dk[1] = sym_qvec[1] - xk[1];
                        dk[2] = sym_qvec[2] - xk[2];
                        dk[0] = dk[0] - std::round(dk[0]);
                        dk[1] = dk[1] - std::round(dk[1]);
                        dk[2] = dk[2] - std::round(dk[2]);
                        if( std::abs(dk[0]) + std::abs(dk[1]) + std::abs(dk[2]) < 1.0e-10 )
                        {
                            find_eq = true;
                            break;
                        }

                        if(Rmg_Symm->time_reversal)
                        {
                            dk[0] = sym_qvec[0] + xk[0];
                            dk[1] = sym_qvec[1] + xk[1];
                            dk[2] = sym_qvec[2] + xk[2];
                            dk[0] = dk[0] - std::round(dk[0]);
                            dk[1] = dk[1] - std::round(dk[1]);
                            dk[2] = dk[2] - std::round(dk[2]);
                            if( std::abs(dk[0]) + std::abs(dk[1]) + std::abs(dk[2]) < 1.0e-10 )
                            {
                                find_eq = true;
                                break;
                            }
                        }
                    }

                    if(find_eq) 
                    {
                        kequiv[i * kmesh[1] * kmesh[2] + j * kmesh[2] + k] = ik;
                        break;
                    }
                }

                if(!find_eq)
                {
                    printf("\n cannot find k equiv  %f %f %f", xk[0], xk[1], xk[2]);
                    throw RmgFatalException() << "kpoint problem" << " in " << __FILE__ << " at line " << __LINE__ << ".\n";
                }

            }

        }
    }

    //  Take the shortest diagonal line as the "shaft" of tetrahedral devision

    double bvec0[3], bvec1[3], bvec2[3];
    for(int i = 0; i < 3; i++)
    {
        bvec0[i] = L.b0[i]/kmesh[0];
        bvec1[i] = L.b1[i]/kmesh[1];
        bvec2[i] = L.b2[i]/kmesh[2];
    }

    std::vector<double> diag_length={0.0, 0.0, 0.0, 0.0};
    double bvec3[3][4];
    for(int i = 0; i < 3; i++)
    {
        bvec3[i][0] = -bvec0[i] + bvec1[i] + bvec2[i];
        bvec3[i][1] =  bvec0[i] - bvec1[i] + bvec2[i];
        bvec3[i][2] =  bvec0[i] + bvec1[i] - bvec2[i];
        bvec3[i][3] =  bvec0[i] + bvec1[i] + bvec2[i];

        diag_length[0] += bvec3[i][0] * bvec3[i][0];
        diag_length[1] += bvec3[i][1] * bvec3[i][1];
        diag_length[2] += bvec3[i][2] * bvec3[i][2];
        diag_length[3] += bvec3[i][3] * bvec3[i][3];

    }

    int id_short = std::min_element(diag_length.begin(), diag_length.end()) - diag_length.begin();

    int ivvec[3][20][6], divvec[4][4], ivvec0[4]={0,0,0,0};
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            divvec[i][j] = 0;
            if(i == j) divvec[i][j] = 1;
        }
    }

    ivvec0[id_short] = 1;
    divvec[id_short][id_short] = -1;

    // set the four corner point of 6 tetra for each cell

    int itet = 0;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if( j == i) continue;
            for(int k = 0; k < 3; k++)
            {
                if(k == i || k == j) continue;
                for(int id = 0; id < 3; id++)
                {
                    ivvec[id][0][itet] = ivvec0[id]; 
                    ivvec[id][1][itet] = ivvec[id][0][itet] + divvec[id][i];
                    ivvec[id][2][itet] = ivvec[id][1][itet] + divvec[id][j];
                    ivvec[id][3][itet] = ivvec[id][2][itet] + divvec[id][k];
                }
                itet++;
            }
        }
    }

    // additional 16 points for each tetra
    for(int id = 0; id < 3; id++)
    {
        for(int it = 0; it < 6; it++)
        {
            ivvec[id][4][it] = 2* ivvec[id][0][it] - ivvec[id][1][it];
            ivvec[id][5][it] = 2* ivvec[id][1][it] - ivvec[id][2][it];
            ivvec[id][6][it] = 2* ivvec[id][2][it] - ivvec[id][3][it];
            ivvec[id][7][it] = 2* ivvec[id][3][it] - ivvec[id][0][it];

            ivvec[id][ 8][it] = 2* ivvec[id][0][it] - ivvec[id][2][it];
            ivvec[id][ 9][it] = 2* ivvec[id][1][it] - ivvec[id][3][it];
            ivvec[id][10][it] = 2* ivvec[id][2][it] - ivvec[id][0][it];
            ivvec[id][11][it] = 2* ivvec[id][3][it] - ivvec[id][1][it];

            ivvec[id][12][it] = 2* ivvec[id][0][it] - ivvec[id][3][it];
            ivvec[id][13][it] = 2* ivvec[id][1][it] - ivvec[id][0][it];
            ivvec[id][14][it] = 2* ivvec[id][2][it] - ivvec[id][1][it];
            ivvec[id][15][it] = 2* ivvec[id][3][it] - ivvec[id][2][it];

            ivvec[id][16][it] = ivvec[id][3][it] - ivvec[id][0][it] + ivvec[id][1][it];
            ivvec[id][17][it] = ivvec[id][0][it] - ivvec[id][1][it] + ivvec[id][2][it];
            ivvec[id][18][it] = ivvec[id][1][it] - ivvec[id][2][it] + ivvec[id][3][it];
            ivvec[id][19][it] = ivvec[id][2][it] - ivvec[id][3][it] + ivvec[id][0][it];
        }
    }

    int itettot = 0;
    for(int i = 0; i < kmesh[0]; i++)
    {
        for(int j = 0; j < kmesh[1]; j++)
        {
            for(int k = 0; k < kmesh[2]; k++)
            {
                for(int it = 0; it < 6; it++)
                {
                    for(int ip = 0; ip < nk_per_tetra; ip++)
                    {
                        int ik0 = i + ivvec[0][ip][itet];
                        int ik1 = j + ivvec[1][ip][itet];
                        int ik2 = k + ivvec[2][ip][itet];
                        ik0 = (kmesh[0] + ik0 % kmesh[0]) %kmesh[0];
                        ik1 = (kmesh[1] + ik1 % kmesh[1]) %kmesh[1];
                        ik2 = (kmesh[2] + ik2 % kmesh[2]) %kmesh[2];
                        int ik = ik0 * kmesh[1] * kmesh[2] + ik1 * kmesh[2] + ik2;
                        tetra[ip][itettot] = kequiv[ik];
                    }

                    itettot ++;
                }
            }
        }
    }
}

void Dos::set_Pij()
{
    int P1[4][4] = { 
        {1440,    0,   30,    0},
        {   0, 1440,    0,   30},
        {  30,    0, 1440,    0},
        {   0,   30,    0, 1440} };
    int P2[4][4] = { 
        { -38,    7,   17,  -28},
        { -28,  -38,    7,   17},
        {  17,  -28,  -38,    7},
        {   7,   17,  -28,  -38} };
    int P3[4][4] = { 
        { -56,    9,  -46,    9},
        {   9,  -56,    9,  -46},
        { -46,    9,  -56,    9},
        {   9,  -46,    9,  -56} };
    int P4[4][4] = { 
        { -38,  -28,   17,    7},
        {   7,  -38,  -28,   17},
        {  17,    7,  -38,  -28},
        { -28,   17,    7,  -38} };
    int P5[4][4] = { 
        { -18,  -18,   12,  -18},
        { -18,  -18,  -18,   12},
        {  12,  -18,  -18,  -18},
        { -18,   12,  -18,  -18} };

    for (int i= 0; i< 4; i++)
        for (int j= 0; j< 4; j++)
        {
            Pij[i][j] = (double)P1[i][j]/1260.0;
            Pij[i][j+4] = (double)P2[i][j]/1260.0;
            Pij[i][j+8] = (double)P3[i][j]/1260.0;
            Pij[i][j+12] = (double)P4[i][j]/1260.0;
            Pij[i][j+16] = (double)P5[i][j]/1260.0;
        }

//    for (int i= 0; i< 4; i++)
//        for (int j= 0; j< 4; j++)
//        {
//            Pij[i][j] = 0.0;
//            if(i == j) Pij[i][j] = 1.0;
//            Pij[i][j+4] = 0.0;
//            Pij[i][j+8] = 0.0;
//            Pij[i][j+12] = 0.0;
//            Pij[i][j+16] = 0.0;
//        }

}
