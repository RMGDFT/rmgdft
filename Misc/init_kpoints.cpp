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

extern "C" int spg_get_ir_reciprocal_mesh(int *grid_address,
                               int map[],
                               const int mesh[3],
                               const int is_shift[3],
                               const int is_time_reversal,
                               const double lattice[3][3],
                               const double *position,
                               const int types[],
                               const int num_atom,
                               const double symprec,
                               const double angprec);

int init_kpoints (int *kmesh, int *kshift)
{

    double *tau;
    const int is_time_reversal = 1;
    int meshsize, num_kpts, count;

    double symprec = 1.0e-5, angprec = 1.0;

    int kpt;
    if(ct.is_gamma)
    {
        num_kpts = 1;
        ct.kp.resize(num_kpts);
        ct.num_kpts = num_kpts;
        ct.kp[0].kpt[ 0 ] = 0.0;
        ct.kp[0].kpt[ 1 ] = 0.0;
        ct.kp[0].kpt[ 2 ] = 0.0;
        ct.kp[0].kweight = 1.0;
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

                bool find_eq (false);
                for(int ik = 0; ik < (int)ct.kp.size(); ik++)
                {
                    for(int isym = 0; isym < Rmg_Symm->nsym; isym++)
                    {
                        sym_qvec[0] = Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                                      Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 0 ] * ct.kp[ik].kpt[1] +
                                      Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 0 ] * ct.kp[ik].kpt[2];
                        sym_qvec[1] = Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 1 ] * ct.kp[ik].kpt[0] +
                                      Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                                      Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 1 ] * ct.kp[ik].kpt[2];
                        sym_qvec[2] = Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 2 ] * ct.kp[ik].kpt[0] +
                                      Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 2 ] * ct.kp[ik].kpt[1] +
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
                }

            }

        }
    }

    ct.num_kpts = ct.kp.size();

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

    }

    if (ct.verbose)
    {
        printf("\n num_k %d", count);
        for(kpt = 0; kpt < num_kpts; kpt++)
            printf("\n kvec %d  %f %f %f %f\n", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);
    }

    return ct.num_kpts;


}                               /* end init_kpoints */

