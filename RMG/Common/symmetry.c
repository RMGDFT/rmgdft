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

#if !GAMMA_PT

static int s[MAX_SYMMETRY][3][3];
static int irg[MAX_SYMMETRY], irt[MAX_IONS][MAX_SYMMETRY];
static int ftau[MAX_SYMMETRY][3], ityp[MAX_IONS];
static rmg_double_t tau[MAX_IONS][3];
static int nsym;



typedef struct
{
    rmg_double_t s[FNZ_GRID][FNY_GRID][FNX_GRID];
} DENS_ARRAY;


/* This routine is used to initialize the symmetry structures */
void init_sym (void)
{
    int nr1, nr2, nr3;
    int ion, kpt, wflag;
    rmg_double_t celldm[6];
    rmg_double_t *xk, *wk;

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


    my_malloc (xk, 3 * ct.num_kpts, rmg_double_t);
    my_malloc (wk, ct.num_kpts, rmg_double_t);

    /* Set up special k-point positions and weights for fortran routines */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        xk[3 * kpt] = ct.kp[kpt].kpt[0];
        xk[3 * kpt + 1] = ct.kp[kpt].kpt[1];
        xk[3 * kpt + 2] = ct.kp[kpt].kpt[2];
        wk[kpt] = ct.kp[kpt].kweight;
    }


    celldm[0] = ct.celldm[0];
    celldm[1] = ct.celldm[1];
    celldm[2] = ct.celldm[2];
    celldm[3] = 0.0;
    celldm[4] = 0.0;
    celldm[5] = 0.0;

    if (ct.boundaryflag != CLUSTER)
    {
        /* Call the symmetry routine */
        symmetry (&ct.ibrav, &s[0][0][0], &nsym, irg, &irt[0][0],
                  &ftau[0][0], &ct.num_ions, &tau[0][0], ityp,
                  &ct.num_kpts, xk, wk, ct.celldm, &nr1, &nr2, &nr3, &wflag);
    }

    my_free (xk);
    my_free (wk);
}                               /* end init_sym */




/* Symmetrizes the density */
void symmetrize_rho (FP0_GRID * rho)
{

    int idx, ix, iy, iz, xoff, yoff, zoff, nr1, nr2, nr3;
    DENS_ARRAY *da;
    rmg_double_t t1;


    /* Wait until all processors arrive at this point */
    /*my_barrier(); */



    /* Get some memory */
    my_malloc (da, 1, DENS_ARRAY);


    for (ix = 0; ix < get_FNX_GRID(); ix++)
        for (iy = 0; iy < get_FNY_GRID(); iy++)
            for (iz = 0; iz < get_FNZ_GRID(); iz++)
                da->s[iz][iy][ix] = 0.0;




    /* Put this processors charge in the correct place */
    pe2xyz (pct.gridpe, &ix, &iy, &iz);
    xoff = ix * pct.FPX0_GRID;
    yoff = iy * pct.FPY0_GRID;
    zoff = iz * pct.FPZ0_GRID;



    for (ix = 0; ix < pct.FPX0_GRID; ix++)
        for (iy = 0; iy < pct.FPY0_GRID; iy++)
            for (iz = 0; iz < pct.FPZ0_GRID; iz++)
                da->s[iz + zoff][iy + yoff][ix + xoff] = rho->s1.b[ix][iy][iz];


    /* Call global sums to give everyone the full array */
    idx = get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID();
    global_sums (&da->s[0][0][0], &idx, pct.grid_comm);


    /* Do the symmetrization on this processor */
    nr1 = get_FNX_GRID();
    nr2 = get_FNY_GRID();
    nr3 = get_FNZ_GRID();

    symrho (&da->s[0][0][0], &nr1, &nr2, &nr3, &nsym, &s[0][0][0], irg, &ftau[0][0]);


    /* Pack density back into correct place */
    t1 = (rmg_double_t) nsym;
    t1 = 1.0 / t1;

    for (ix = 0; ix < pct.FPX0_GRID; ix++)
        for (iy = 0; iy < pct.FPY0_GRID; iy++)
            for (iz = 0; iz < pct.FPZ0_GRID; iz++)
                rho->s1.b[ix][iy][iz] = da->s[iz + zoff][iy + yoff][ix + xoff] * t1;


    /* Release our memory */
    my_free (da);


}                               /* end symmetrize_rho */



/* This routine is used to symmetrize the forces (MBN)*/
void symforce (void)
{
    int ion, nr1, nr2, nr3;
    rmg_double_t celldm[6], force[MAX_IONS][3];

    nr1 = get_NX_GRID();
    nr2 = get_NY_GRID();
    nr3 = get_NZ_GRID();

    celldm[0] = ct.celldm[0];
    celldm[1] = ct.celldm[1];
    celldm[2] = ct.celldm[2];
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

    fsymforces (&force[0][0], &s[0][0][0], irg, &irt[0][0], &ct.num_ions, &ct.ibrav,
                &nsym, celldm, &nr1, &nr2, &nr3);

    /* Store forces back in c format */

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        ct.ions[ion].force[ct.fpt[0]][0] = force[ion][0];
        ct.ions[ion].force[ct.fpt[0]][1] = force[ion][1];
        ct.ions[ion].force[ct.fpt[0]][2] = force[ion][2];
    }

}                               /* end symforce */

#endif

/******/
