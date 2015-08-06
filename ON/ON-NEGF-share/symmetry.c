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
#include "prototypes_on.h"

#if !GAMMA_PT

static int s[MAX_SYMMETRY][3][3];
static int irg[MAX_SYMMETRY], irt[MAX_IONS][MAX_SYMMETRY];
static int ftau[MAX_SYMMETRY][3], ityp[MAX_IONS];
static double tau[MAX_IONS][3], xk[MAX_KPTS][3], wk[MAX_KPTS];
static int nsym;



/* This routine is used to initialize the symmetry structures */
void init_sym (void)
{
    int nr1, nr2, nr3;
    int ion, kpt, wflag;

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


    /* Set up special k-point positions and weights for fortran routines */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        xk[kpt][0] = ct.kp[kpt].kpt[0];
        xk[kpt][1] = ct.kp[kpt].kpt[1];
        xk[kpt][2] = ct.kp[kpt].kpt[2];
        wk[kpt] = ct.kp[kpt].kweight;
    }


    if (ct.boundaryflag != CLUSTER)
    {
        int ibrav = get_ibrav_type();
        /* Call the symmetry routine */
        symmetry (&ibrav, &s[0][0][0], &nsym, irg, &irt[0][0],
                  &ftau[0][0], &ct.num_ions, &tau[0][0], ityp,
                  &ct.num_kpts, &xk[0][0], wk, ct.celldm, &nr1, &nr2, &nr3, &wflag);
    }

}                               /* end init_sym */




/* Symmetrizes the density */
void symmetrize_rho (double * rho)
{

    int idx, ix, iy, iz, xoff, yoff, zoff, nr1, nr2, nr3;
    double *da;
    double t1;
    int idx1, incx, incy;
    /* Wait until all processors arrive at this point */
    my_barrier ();



    /* Get some memory */
    my_malloc_init( da, ct.psi_nbasis, double );

    idx = 0;
    for (ix = 0; ix < get_NX_GRID(); ix++)
    {
        for (iy = 0; iy < get_NY_GRID(); iy++)
        {
            for (iz = 0; iz < get_NZ_GRID(); iz++)
            {
                da[idx] = 0.0;
                idx++;
            }
        }
    }



    /* Put this processors charge in the correct place */
    pe2xyz (pct.gridpe, &ix, &iy, &iz);
    xoff = ix * get_PX0_GRID();
    yoff = iy * get_PY0_GRID();
    zoff = iz * get_PZ0_GRID();


    incx = get_NY_GRID() * get_NZ_GRID();
    incy = get_NZ_GRID();

    idx = 0;
    for (ix = 0; ix < get_PX0_GRID(); ix++)
    {
        for (iy = 0; iy < get_PY0_GRID(); iy++)
        {
            for (iz = 0; iz < get_PZ0_GRID(); iz++)
            {
                idx1 = (ix + xoff) * incx + (iy + yoff) * incy + (iz + zoff);

                da[idx1] = rho[idx];
                idx++;
            }
        }
    }

    /* Call global sums to give everyone the full array */
    idx = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
    global_sums (da, &idx, pct.grid_comm);


    /* Do the symmetrization on this processor */
    nr1 = get_NX_GRID();
    nr2 = get_NY_GRID();
    nr3 = get_NZ_GRID();

    symrho (da, &nr1, &nr2, &nr3, &nsym, &s[0][0][0], irg, &ftau[0][0]);


    /* Pack density back into correct place */
    t1 = (double) nsym;
    t1 = 1.0 / t1;

    idx = 0;
    for (ix = 0; ix < get_PX0_GRID(); ix++)
    {
        for (iy = 0; iy < get_PY0_GRID(); iy++)
        {
            for (iz = 0; iz < get_PZ0_GRID(); iz++)
            {
                idx1 = (ix + xoff) * incx + (iy + yoff) * incy + (iz + zoff);

                rho[idx] = da[idx1] * t1;

                idx++;
            }
        }
    }


    /* Release our memory */
    my_free(da);


}                               /* end symmetrize_rho */



/* This routine is used to symmetrize the forces (MBN)*/
void symforce (void)
{
    int ion, nr1, nr2, nr3;
    double celldm[6], force[MAX_IONS][3];

    nr1 = get_NX_GRID();
    nr2 = get_NY_GRID();
    nr3 = get_NZ_GRID();

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

    int ibrav = get_ibrav_type();
    fsymforces (&force[0][0], &s[0][0][0], irg, &irt[0][0], &ct.num_ions, &ibrav,
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
