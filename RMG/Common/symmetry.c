/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/symmetry.c *****
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
 *   C wrapper routines to call fortran functions which actually do the
 *   symmetrization.
 * 
 *   void init_sym(void)
 *     This routine is used to initialize the symmetry structures 
 *   void symmetrize_rho(P0_GRID *rho)
 *     Symmetrizes the density
 *   void symforce(void)
 *     This routine is used to symmetrize the forces (MBN)
 * INPUTS
 *   nothing
 * OUTPUT
 *   charge density and forces are updated 
 * PARENTS
 *   force.c get_rho.c init.c
 * CHILDREN
 *   to_crystal.c symmetry.f
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include "main.h"


static int s[MAX_SYMMETRY][3][3];
static int irg[MAX_SYMMETRY], irt[MAX_IONS][MAX_SYMMETRY];
static int ftau[MAX_SYMMETRY][3], ityp[MAX_IONS];
static double tau[MAX_IONS][3];
static int nsym;
static double translation[MAX_SYMMETRY][3];



/* This routine is used to initialize the symmetry structures */
void init_sym (void)
{
    int nr1, nr2, nr3;
    int ion, kpt, wflag;
    double celldm[6];
    double *xk, *wk;
    int ibrav = get_ibrav_type();
    double a0[3], a1[3], a2[3], omega;

    double lattice[3][3];

    lattice[0][0] = get_a0(0);
    lattice[0][1] = get_a0(1);
    lattice[0][2] = get_a0(2);
    lattice[1][0] = get_a1(0);
    lattice[1][1] = get_a1(1);
    lattice[1][2] = get_a1(2);
    lattice[2][0] = get_a2(0);
    lattice[2][1] = get_a2(1);
    lattice[2][2] = get_a2(2);

    a0[0] = get_a0(0);
    a0[1] = get_a0(1);
    a0[2] = get_a0(2);
    a1[0] = get_a1(0);
    a1[1] = get_a1(1);
    a1[2] = get_a1(2);
    a2[0] = get_a2(0);
    a2[1] = get_a2(1);
    a2[2] = get_a2(2);
    omega = get_omega();

    /* This function uses MAX_IONS as a limit for array sizes.
     * It is, of course, possible to allocate these arrays dynamically,
     * but I am not sure if this would not break FORTRAN routines to which
     * those arrays are passed*/

    if (ct.num_ions >= MAX_IONS)
        error_handler ("Too many ions, increase MAX_IONS in params.h");

    nr1 = get_NX_GRID();
    nr2 = get_NY_GRID();
    nr3 = get_NZ_GRID();

    /* Only have PE zero output symmetry information */
    wflag = pct.gridpe;


    /* Set up atomic positions and species for fortran routines */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        to_crystal (&tau[ion][0], ct.ions[ion].crds);
        ityp[ion] = ct.ions[ion].species;
    }


    my_malloc (xk, 3 * ct.num_kpts, double);
    my_malloc (wk, ct.num_kpts, double);

    /* Set up special k-point positions and weights for fortran routines */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        xk[3 * kpt] = ct.kp[kpt].kpt[0];
        xk[3 * kpt + 1] = ct.kp[kpt].kpt[1];
        xk[3 * kpt + 2] = ct.kp[kpt].kpt[2];
        wk[kpt] = ct.kp[kpt].kweight;
    }


    celldm[0] = get_celldm(0);
    celldm[1] = get_celldm(1);
    celldm[2] = get_celldm(2);
    celldm[3] = 0.0;
    celldm[4] = 0.0;
    celldm[5] = 0.0;

    nsym = spg_get_symmetry(&s[0][0][0], translation,  MAX_SYMMETRY, lattice, tau, ityp, ct.num_ions, 1e-5);
    printf("\n nysm %d", nsym);
    int i, j;
    for(kpt = 0; kpt < nsym; kpt++)
    {
        printf("\n sym operation # %d", kpt);
            for(i = 0; i < 3; i++)
            {
                printf("\n");
                for(j = 0; j < 3; j++)
                {
                    printf(" %d ", s[kpt][i][j]);
                }
            }
        printf("\n %f %f %f ", translation[kpt][0], translation[kpt][1], translation[kpt][2]);
    }

    if (ct.boundaryflag != CLUSTER)
    {
        /* Call the symmetry routine */
        symmetry (&ibrav, &s[0][0][0], &nsym, irg, &irt[0][0],
                &ftau[0][0], &ct.num_ions, &tau[0][0], ityp,
                &ct.num_kpts, xk, wk, celldm, &nr1, &nr2, &nr3, a0, a1, a2, &omega, &wflag);
    }

    my_free (xk);
    my_free (wk);
}                               /* end init_sym */




/* Symmetrizes the density */
void symmetrize_rho (double * rho)
{

    int idx, ix, iy, iz, xoff, yoff, zoff;
    double *da;
    double t1;

    int FPX0_GRID = get_FPX0_GRID();
    int FPY0_GRID = get_FPY0_GRID();
    int FPZ0_GRID = get_FPZ0_GRID();
    int FP0_BASIS = get_FP0_BASIS();
    int FNX_GRID = get_FNX_GRID();
    int FNY_GRID = get_FNY_GRID();
    int FNZ_GRID = get_FNZ_GRID();
    int FN_BASIS = FNX_GRID * FNY_GRID * FNZ_GRID;

    /* Wait until all processors arrive at this point */
    /*my_barrier(); */



    /* Get some memory */
    my_malloc (da, FN_BASIS, double);


    for(idx = 0;idx < FN_BASIS;idx++) {
        da[idx] = 0.0;
    }



    /* Put this processors charge in the correct place */
    pe2xyz (pct.gridpe, &ix, &iy, &iz);
    xoff = ix * FPX0_GRID;
    yoff = iy * FPY0_GRID;
    zoff = iz * FPZ0_GRID;

    int incx = FPY0_GRID * FPZ0_GRID;
    int incy = FPZ0_GRID;

    int incx1 = 1;
    int incy1 = FNX_GRID;
    int incz1 = FNX_GRID * FNY_GRID;

    for (ix = 0; ix < FPX0_GRID; ix++) {
        for (iy = 0; iy < FPY0_GRID; iy++) {
            for (iz = 0; iz < FPZ0_GRID; iz++) {
                da[(iz + zoff)*incz1 + (iy + yoff)*incy1 + (ix + xoff)*incx1] = rho[ix * incx + iy*incy + iz];
            }
        }
    }

    /* Call global sums to give everyone the full array */
    global_sums (da, &FN_BASIS, pct.grid_comm);


    /* Do the symmetrization on this processor */

    symrho (da, &FNX_GRID, &FNY_GRID, &FNZ_GRID, &nsym, &s[0][0][0], irg, &ftau[0][0]);


    /* Pack density back into correct place */
    t1 = (double) nsym;
    t1 = 1.0 / t1;

    for (ix = 0; ix < FPX0_GRID; ix++) {
        for (iy = 0; iy < FPY0_GRID; iy++) {
            for (iz = 0; iz < FPZ0_GRID; iz++) {
                rho[ix * incx + iy*incy + iz] = da[(iz + zoff)*incz1 + (iy + yoff)*incy1 + (ix + xoff)*incx1] * t1;
            }
        }
    }


    /* Release our memory */
    my_free (da);


}                               /* end symmetrize_rho */



/* This routine is used to symmetrize the forces (MBN)*/
void symforce (void)
{
    int ion;
    double celldm[6], force[MAX_IONS][3];

    int NX_GRID = get_NX_GRID();
    int NY_GRID = get_NY_GRID();
    int NZ_GRID = get_NZ_GRID();

    int ibrav = get_ibrav_type();

    celldm[0] = get_celldm(0);
    celldm[1] = get_celldm(1);
    celldm[2] = get_celldm(2);
    celldm[3] = 0.0;
    celldm[4] = 0.0;
    celldm[5] = 0.0;

    /* Convert forces to format expected by fortran routines */

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        force[ion][0] = ct.ions[ion].force[ct.fpt[0]][0];
        force[ion][1] = ct.ions[ion].force[ct.fpt[0]][1];
        force[ion][2] = ct.ions[ion].force[ct.fpt[0]][2];
    }

    fsymforces (&force[0][0], &s[0][0][0], irg, &irt[0][0], &ct.num_ions, &ibrav,
            &nsym, celldm, &NX_GRID, &NY_GRID, &NZ_GRID);

    /* Store forces back in c format */

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        ct.ions[ion].force[ct.fpt[0]][0] = force[ion][0];
        ct.ions[ion].force[ct.fpt[0]][1] = force[ion][1];
        ct.ions[ion].force[ct.fpt[0]][2] = force[ion][2];
    }

}                               /* end symforce */

/******/
