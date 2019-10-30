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

int init_kpoints (int *mesh, int *is_shift)
{

    double *tau;
    int *grid_address, *map, *weight;
    const int is_time_reversal = 1;
    int meshsize, num_kpts, count;

    int ion, *ityp;
    double symprec = 1.0e-5, angprec = 1.0;

    int kpt;
    if(mesh[0] == 1 && mesh[1] == 1 && mesh[2] == 1 
            && is_shift[0] == 0 && is_shift[1] == 0 && is_shift[2] == 0)
    {
        num_kpts = 1;
        ct.kp.resize(num_kpts);
        ct.num_kpts = num_kpts;
        ct.kp[0].kpt[ 0 ] = 0.0;
        ct.kp[0].kpt[ 1 ] = 0.0;
        ct.kp[0].kpt[ 2 ] = 0.0;
        ct.kp[0].kweight = 1.0;
        ct.is_gamma = 1;
        ct.is_use_symmetry = false;
        return ct.num_kpts;
    }

    double lattice[3][3];
    lattice[0][0] = get_a0(0);
    lattice[1][0] = get_a0(1);
    lattice[2][0] = get_a0(2);
    lattice[0][1] = get_a1(0);
    lattice[1][1] = get_a1(1);
    lattice[2][1] = get_a1(2);
    lattice[0][2] = get_a2(0);
    lattice[1][2] = get_a2(1);
    lattice[2][2] = get_a2(2);


    meshsize = mesh[0] * mesh[1] * mesh[2];
    grid_address = new int[meshsize*5];

    map = grid_address + 3*meshsize;
    weight = map + meshsize;

    tau = new double[ct.num_ions*3];
    ityp = new int[ct.num_ions];

    /* Set up atomic positions and species for fortran routines */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        to_crystal (&tau[ion*3], Atoms[ion].crds);
        ityp[ion] = Atoms[ion].species;
    }





    num_kpts = spg_get_ir_reciprocal_mesh(grid_address, map,mesh, is_shift, 
            is_time_reversal, lattice, tau, ityp, ct.num_ions, symprec, angprec);

    //    printf("\n num_k %d", num_kpts);
    //    for(kpt = 0; kpt < meshsize; kpt++)
    //        printf("\n kvec %d  %d %d %d %d", kpt, grid_address[kpt*3 +0], grid_address[kpt*3 +1], grid_address[kpt*3 +2], map[kpt]);


    ct.kp.resize(num_kpts);
    ct.num_kpts = num_kpts;
    for(kpt = 0; kpt < meshsize; kpt++)
        weight[kpt] = 0;
    for(kpt = 0; kpt < meshsize; kpt++)
        weight[map[kpt]]++;

    count = 0;
    for(kpt = 0; kpt < meshsize; kpt++)
    {
        if(weight[kpt] != 0)
        {
            ct.kp[count].kpt[ 0 ] = (grid_address[kpt*3 + 0] + is_shift[0]*0.5)/mesh[0];
            ct.kp[count].kpt[ 1 ] = (grid_address[kpt*3 + 1] + is_shift[1]*0.5)/mesh[1];
            ct.kp[count].kpt[ 2 ] = (grid_address[kpt*3 + 2] + is_shift[2]*0.5)/mesh[2];
            ct.kp[count].kweight = (double)weight[kpt]/meshsize;
            count++;
        }
    }

    assert(count == num_kpts);
    /*Not necessary and it ends up being the first thing printed to stdout*/

    ct.is_gamma = true;
    for (int kpt = 0; kpt < ct.num_kpts; kpt++) {
        double v1, v2, v3;
        v1 = 0.0;
        v2 = 0.0;
        v3 = 0.0;

        for(int ir = 0; ir<3; ir++)
        {
            v1 = ct.kp[kpt].kpt[0] *Rmg_L.b0[0]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[0] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[0];
            v2 = ct.kp[kpt].kpt[0] *Rmg_L.b0[1]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[1] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[1];
            v3 = ct.kp[kpt].kpt[0] *Rmg_L.b0[2]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[2] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[2];
        }

        ct.kp[kpt].kvec[0] = v1 * twoPI;
        ct.kp[kpt].kvec[1] = v2 * twoPI;
        ct.kp[kpt].kvec[2] = v3 * twoPI;
        ct.kp[kpt].kmag = (v1 * v1 + v2 * v2 + v3 * v3) * twoPI * twoPI;

        if(ct.kp[kpt].kmag != 0.0) ct.is_gamma = false;
    }

    if(ct.is_gamma) 
    {
        ct.is_use_symmetry = 0;
    }

    if (ct.verbose)
    {
        printf("\n num_k %d", count);
        for(kpt = 0; kpt < num_kpts; kpt++)
          printf("\n kvec %d  %f %f %f %f\n", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);
    }

    delete [] grid_address;
    delete [] tau;
    delete [] ityp;

    return num_kpts;


}                               /* end init_kpoints */

