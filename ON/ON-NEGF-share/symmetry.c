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
#include "main_on.h"



static int s[MAX_SYMMETRY][3][3];
static int irg[MAX_SYMMETRY], irt[MAX_IONS][MAX_SYMMETRY];
static int ftau[MAX_SYMMETRY][3], ityp[MAX_IONS];
static rmg_double_t tau[MAX_IONS][3], xk[MAX_KPTS][3], wk[MAX_KPTS];
static int nsym;



/* This routine is used to initialize the symmetry structures */
void init_sym (void)
{
    int nr1, nr2, nr3;
    int ion, kpt, wflag;

    nr1 = ct.psi_nxgrid;
    nr2 = ct.psi_nygrid;
    nr3 = ct.psi_nzgrid;

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
        /* Call the symmetry routine */
        symmetry (&ct.ibrav, &s[0][0][0], &nsym, irg, &irt[0][0],
                  &ftau[0][0], &ct.num_ions, &tau[0][0], ityp,
                  &ct.num_kpts, &xk[0][0], wk, ct.celldm, &nr1, &nr2, &nr3, &wflag);
    }

}                               /* end init_sym */




/* Symmetrizes the density */
void symmetrize_rho (rmg_double_t * rho)
{

    int idx, ix, iy, iz, xoff, yoff, zoff, nr1, nr2, nr3;
    rmg_double_t *da;
    rmg_double_t t1;
    int idx1, incx, incy;
    /* Wait until all processors arrive at this point */
    my_barrier ();



    /* Get some memory */
    my_malloc_init( da, ct.psi_nbasis, rmg_double_t );

    idx = 0;
    for (ix = 0; ix < ct.psi_nxgrid; ix++)
    {
        for (iy = 0; iy < ct.psi_nygrid; iy++)
        {
            for (iz = 0; iz < ct.psi_nzgrid; iz++)
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
    idx = ct.psi_nxgrid * ct.psi_nygrid * ct.psi_nzgrid;
    global_sums (da, &idx, pct.grid_comm);


    /* Do the symmetrization on this processor */
    nr1 = ct.psi_nxgrid;
    nr2 = ct.psi_nygrid;
    nr3 = ct.psi_nzgrid;

    symrho (da, &nr1, &nr2, &nr3, &nsym, &s[0][0][0], irg, &ftau[0][0]);


    /* Pack density back into correct place */
    t1 = (rmg_double_t) nsym;
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
    rmg_double_t celldm[6], force[MAX_IONS][3];

    nr1 = ct.psi_nxgrid;
    nr2 = ct.psi_nygrid;
    nr3 = ct.psi_nzgrid;

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

/******/
