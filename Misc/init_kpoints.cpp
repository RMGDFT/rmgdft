/****f* QMD-MGDFT/init_pos.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void init_pos()
 *   if positions are read in cartesian coordinates, get crystal coordinates
 *   if positions are read in crystal coordinates, get cartesian coordinates
 * INPUTS
 *   nothing
 * OUTPUT
 *   coordinates are stored in ct.ion 
 * PARENTS
 *   init.c
 * CHILDREN
 *   to_cartesian.c
 * SOURCE
 */

#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "transition.h"
#include "rmg_mangling.h"
#include "qk_compat.h"

extern "C" int *s_fortran_ptr;

int init_kpoints (int *kmesh, int *kshift)
{
    int max_kpts = kmesh[0] * kmesh[1] * kmesh[2];
    int magnetic_sym = 0;
    int no_z_inv = 0;


    //double *tau;
    int num_kpts=1;

    int kpt;
    for(int i = 0; i < 3; i++)
    {
        ct.klist.kpoint_mesh[i] = kmesh[i];
        ct.klist.kpoint_is_shift[i] = kshift[i];
    }

    ct.klist.num_k_all = kmesh[0] * kmesh[1] * kmesh[2];
    ct.klist.k_all_xtal.resize(boost::extents[ct.klist.num_k_all][3]);
    ct.klist.k_all_cart.resize(boost::extents[ct.klist.num_k_all][3]);
    ct.klist.k_map_index.resize(ct.klist.num_k_all, 0);
    ct.klist.k_map_symm.resize(ct.klist.num_k_all, 0);

    ct.klist.num_k_ext = (kmesh[0] + 2) * (kmesh[1]+2) * (kmesh[2]+2) - ct.klist.num_k_all;
    ct.klist.k_ext_xtal.resize(boost::extents[ct.klist.num_k_ext][3]);
    ct.klist.k_ext_cart.resize(boost::extents[ct.klist.num_k_ext][3]);
    
    if(ct.is_gamma)
    {
        num_kpts = 1;
        ct.kp.resize(num_kpts);
        ct.num_kpts = num_kpts;
        ct.kp[0].kpt[ 0 ] = 0.0;
        ct.kp[0].kpt[ 1 ] = 0.0;
        ct.kp[0].kpt[ 2 ] = 0.0;
        ct.kp[0].kweight = 1.0;

        ct.klist.num_k_ire = 1;
        ct.klist.k_ire_xtal.resize(boost::extents[ct.klist.num_k_ire][3]);
        ct.klist.k_ire_cart.resize(boost::extents[ct.klist.num_k_ire][3]);

        for(int i = 0; i < 3; i++)
        {
            ct.klist.k_ire_xtal[0][i] = 0.0;
            ct.klist.k_all_xtal[0][i] = 0.0;
            ct.klist.k_ire_cart[0][i] = 0.0;
            ct.klist.k_all_cart[0][i] = 0.0;
        }
    
        ct.klist.kweight.resize(ct.klist.num_k_ire, 1.0);

        return ct.num_kpts;
    }


    ct.kp.clear();

    double weight_one = 1.0/(kmesh[0] * kmesh[1] * kmesh[2]);

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

                if(xk[0] > 0.5) xk[0] -= 1.0;
                if(xk[1] > 0.5) xk[1] -= 1.0;
                if(xk[2] > 0.5) xk[2] -= 1.0;


                int idx = i * kmesh[1] * kmesh[2] + j * kmesh[2] + k;
                ct.klist.k_all_xtal[idx][0] = xk[0];
                ct.klist.k_all_xtal[idx][1] = xk[1];
                ct.klist.k_all_xtal[idx][2] = xk[2];
                bool find_eq (false);
                for(int ik = 0; ik < (int)ct.kp.size(); ik++)
                {
                    int isym;
                    for(isym = 0; isym < Rmg_Symm->nsym_full; isym++)
                    {
                        sym_qvec[0] = Rmg_Symm->full_sym_rotate[isym * 9 + 0 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                                      Rmg_Symm->full_sym_rotate[isym * 9 + 1 * 3 + 0 ] * ct.kp[ik].kpt[1] +
                                      Rmg_Symm->full_sym_rotate[isym * 9 + 2 * 3 + 0 ] * ct.kp[ik].kpt[2];
                        sym_qvec[1] = Rmg_Symm->full_sym_rotate[isym * 9 + 0 * 3 + 1 ] * ct.kp[ik].kpt[0] +
                                      Rmg_Symm->full_sym_rotate[isym * 9 + 1 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                                      Rmg_Symm->full_sym_rotate[isym * 9 + 2 * 3 + 1 ] * ct.kp[ik].kpt[2];
                        sym_qvec[2] = Rmg_Symm->full_sym_rotate[isym * 9 + 0 * 3 + 2 ] * ct.kp[ik].kpt[0] +
                                      Rmg_Symm->full_sym_rotate[isym * 9 + 1 * 3 + 2 ] * ct.kp[ik].kpt[1] +
                                      Rmg_Symm->full_sym_rotate[isym * 9 + 2 * 3 + 2 ] * ct.kp[ik].kpt[2];

                        if(Rmg_Symm->time_rev[isym])
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
                                ct.klist.k_map_index[idx] = ik;
                                ct.klist.k_map_symm[idx] = -(isym+1);
                                ct.klist.k_all_xtal[idx][0] = -sym_qvec[0];
                                ct.klist.k_all_xtal[idx][1] = -sym_qvec[1];
                                ct.klist.k_all_xtal[idx][2] = -sym_qvec[2];
                                break;
                            }

                        }
                        else 
                        {
                            dk[0] = sym_qvec[0] - xk[0];
                            dk[1] = sym_qvec[1] - xk[1];
                            dk[2] = sym_qvec[2] - xk[2];
                            dk[0] = dk[0] - std::round(dk[0]);
                            dk[1] = dk[1] - std::round(dk[1]);
                            dk[2] = dk[2] - std::round(dk[2]);
                            if( std::abs(dk[0]) + std::abs(dk[1]) + std::abs(dk[2]) < 1.0e-10 )
                            {
                                find_eq = true;
                                ct.klist.k_map_index[idx] = ik;
                                ct.klist.k_map_symm[idx] = isym + 1;
                                ct.klist.k_all_xtal[idx][0] =  sym_qvec[0];
                                ct.klist.k_all_xtal[idx][1] =  sym_qvec[1];
                                ct.klist.k_all_xtal[idx][2] =  sym_qvec[2];
                                break;
                            }

                            if(Rmg_Symm->time_reversal  )
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
                                    ct.klist.k_map_index[idx] = ik;
                                    ct.klist.k_map_symm[idx] = -(isym+1);
                                    ct.klist.k_all_xtal[idx][0] = -sym_qvec[0];
                                    ct.klist.k_all_xtal[idx][1] = -sym_qvec[1];
                                    ct.klist.k_all_xtal[idx][2] = -sym_qvec[2];
                                    break;
                                }
                            }
                        }
                    }

                    if(find_eq) 
                    {
                        ct.kp[ik].kweight += weight_one;
                        break;
                    }
                }

                if(!find_eq)
                {
                    KSTRUCT onekp;
                    onekp.kpt[0] = xk[0];
                    onekp.kpt[1] = xk[1];
                    onekp.kpt[2] = xk[2];
                    onekp.kweight = weight_one;
                    ct.kp.push_back(onekp);
                    ct.klist.k_map_index[idx] = ct.kp.size()-1;
                    ct.klist.k_map_symm[idx] = 1; 
                }

            }

        }
    }

///////////////////////////////////////
    int nks = ct.kp.size();
    std::vector<double> xk1;
    xk1.resize(3*max_kpts);
    std::vector<double> wk;
    wk.resize(max_kpts);
    //printf("nnnnnsym = %d  npk = %d  nkstot = %d  nrot = %d\n",Rmg_Symm->nsym, max_kpts, nks, Rmg_Symm->nsym_full);
    // irreducible_bz uses these units
    xk1.clear();
    xk1.resize(3*max_kpts);
    int ioff = 0;
    for (auto it = begin (ct.kp); it != end (ct.kp); ++it) {
        double v1, v2, v3;

        v1 = it->kpt[0] *Rmg_L.b0[0]
            + it->kpt[1] *Rmg_L.b1[0]
            + it->kpt[2] *Rmg_L.b2[0];
        v2 = it->kpt[0] *Rmg_L.b0[1]
            + it->kpt[1] *Rmg_L.b1[1]
            + it->kpt[2] *Rmg_L.b2[1];
        v3 = it->kpt[0] *Rmg_L.b0[2]
            + it->kpt[1] *Rmg_L.b1[2]
            + it->kpt[2] *Rmg_L.b2[2];

        xk1[ioff*3 + 0] = v1*Rmg_L.celldm[0];
        xk1[ioff*3 + 1] = v2*Rmg_L.celldm[0];
        xk1[ioff*3 + 2] = v3*Rmg_L.celldm[0];
        wk[ioff] = it->kweight;
        ioff++;
    }
    if(ct.verbose && pct.gridpe == 0) 
       printf("II1 nrot = %d nsym = %d  nks = %d\n",Rmg_Symm->nsym_full, Rmg_Symm->nsym, nks);

    // irreducible_bz is cartesian
    irreducible_bz( &Rmg_Symm->nsym_full,
                    Rmg_Symm->full_sym_rotate.data(),
                    &Rmg_Symm->nsym,
                    &ct.time_reversal,
                    &magnetic_sym,
                    Rmg_L.at,         // Use Rmg_L.at
                    Rmg_L.bg,         // Use Rmg_L.bg
                    &max_kpts,        // Max number of k-points
                    &nks,             // Number of k-points
                    xk1.data(),         // K-points stored as triplets
                    wk.data(),         // Their weights
                    Rmg_Symm->time_rev.data());       // Use Rmg_Symm.time_rev.data()

//printf("irreducible_bz  NROT = %d  NSYM = %d  NKS = %d  MAGSYM = %d  MINUS_Q = %d\n",
//Rmg_Symm->nsym_full, Rmg_Symm->nsym, nks, magnetic_sym, minus_q);

    ct.kp.clear();
    ct.kp.resize(nks);
    for(int ik=0;ik < nks;ik++)
    {
        ct.kp[ik].kvec[0] = xk1[ik*3 + 0];
        ct.kp[ik].kvec[1] = xk1[ik*3 + 1];
        ct.kp[ik].kvec[2] = xk1[ik*3 + 2];
        ct.kp[ik].kweight = wk[ik];
    double v1, v2, v3;

    v1 =  Rmg_L.a0[0] * xk1[ik*3+0] 
        + Rmg_L.a0[1] * xk1[ik*3+1] 
        + Rmg_L.a0[2] * xk1[ik*3+2];
    v2 =  Rmg_L.a1[0] * xk1[ik*3+0] 
        + Rmg_L.a1[1] * xk1[ik*3+1] 
        + Rmg_L.a1[2] * xk1[ik*3+2];
    v3 =  Rmg_L.a2[0] * xk1[ik*3+0] 
        + Rmg_L.a2[1] * xk1[ik*3+1] 
        + Rmg_L.a2[2] * xk1[ik*3+2];
        ct.kp[ik].kpt[0] = v1/Rmg_L.celldm[0];
        ct.kp[ik].kpt[1] = v2/Rmg_L.celldm[0];
        ct.kp[ik].kpt[2] = v3/Rmg_L.celldm[0];
//printf("JJJJ0  %f  %f  %f\n", xk1[ik*3 + 0], xk1[ik*3 + 1], xk1[ik*3 + 2]);
//printf("JJJJ1  %f  %f  %f\n", v1/Rmg_L.celldm[0], v2/Rmg_L.celldm[0], v3/Rmg_L.celldm[0]);
    }

///////////////////////////////////////

    ct.num_kpts = ct.kp.size();
    ct.klist.num_k_ire = ct.kp.size();
    ct.klist.k_ire_xtal.resize(boost::extents[ct.klist.num_k_ire][3]);
    ct.klist.k_ire_cart.resize(boost::extents[ct.klist.num_k_ire][3]);
    ct.klist.kweight.resize(ct.klist.num_k_ire, 1.0);

    /*Not necessary and it ends up being the first thing printed to stdout*/

    for (int kpt = 0; kpt < ct.num_kpts; kpt++) {
        double v1, v2, v3;

        v1 = ct.kp[kpt].kpt[0] *Rmg_L.b0[0]
            + ct.kp[kpt].kpt[1] *Rmg_L.b1[0] 
            + ct.kp[kpt].kpt[2] *Rmg_L.b2[0];
        v2 = ct.kp[kpt].kpt[0] *Rmg_L.b0[1]
            + ct.kp[kpt].kpt[1] *Rmg_L.b1[1] 
            + ct.kp[kpt].kpt[2] *Rmg_L.b2[1];
        v3 = ct.kp[kpt].kpt[0] *Rmg_L.b0[2]
            + ct.kp[kpt].kpt[1] *Rmg_L.b1[2] 
            + ct.kp[kpt].kpt[2] *Rmg_L.b2[2];

        ct.kp[kpt].kvec[0] = v1 * twoPI;
        ct.kp[kpt].kvec[1] = v2 * twoPI;
        ct.kp[kpt].kvec[2] = v3 * twoPI;
        ct.kp[kpt].kmag = (v1 * v1 + v2 * v2 + v3 * v3) * twoPI * twoPI;

        ct.klist.k_ire_xtal[kpt][0] = ct.kp[kpt].kpt[0];
        ct.klist.k_ire_xtal[kpt][1] = ct.kp[kpt].kpt[1];
        ct.klist.k_ire_xtal[kpt][2] = ct.kp[kpt].kpt[2];
        ct.klist.k_ire_cart[kpt][0] = ct.kp[kpt].kvec[0];
        ct.klist.k_ire_cart[kpt][1] = ct.kp[kpt].kvec[1];
        ct.klist.k_ire_cart[kpt][2] = ct.kp[kpt].kvec[2];
        ct.klist.kweight[kpt] = ct.kp[kpt].kweight;

    }

    for (int kpt = 0; kpt < ct.klist.num_k_all; kpt++) {
        double v1, v2, v3;

        v1 = ct.klist.k_all_xtal[kpt][0] *Rmg_L.b0[0]
            + ct.klist.k_all_xtal[kpt][1] *Rmg_L.b1[0] 
            + ct.klist.k_all_xtal[kpt][2] *Rmg_L.b2[0];
        v2 = ct.klist.k_all_xtal[kpt][0] *Rmg_L.b0[1]
            + ct.klist.k_all_xtal[kpt][1] *Rmg_L.b1[1] 
            + ct.klist.k_all_xtal[kpt][2] *Rmg_L.b2[1];
        v3 = ct.klist.k_all_xtal[kpt][0] *Rmg_L.b0[2]
            + ct.klist.k_all_xtal[kpt][1] *Rmg_L.b1[2] 
            + ct.klist.k_all_xtal[kpt][2] *Rmg_L.b2[2];
        ct.klist.k_all_cart[kpt][0] = v1 * twoPI;
        ct.klist.k_all_cart[kpt][1] = v2 * twoPI;
        ct.klist.k_all_cart[kpt][2] = v3 * twoPI;
    }

    if (ct.verbose)
    {
        printf("\n num_k %d", ct.num_kpts);
        for(kpt = 0; kpt < ct.num_kpts; kpt++)
            printf("\n kvec %d  %f %f %f %f\n", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);
        for(kpt = 0; kpt < ct.klist.num_k_all; kpt++)
            printf("\n kall %d %f %f %f %d %d %d", kpt,
                    ct.klist.k_all_cart[kpt][0],ct.klist.k_all_cart[kpt][1],ct.klist.k_all_cart[kpt][2],ct.klist.k_map_index[kpt],ct.klist.k_map_symm[kpt],
                    (int)Rmg_Symm->time_rev[std::abs(ct.klist.k_map_symm[kpt])-1 ]);
    }

    return ct.num_kpts;


}                               /* end init_kpoints */

