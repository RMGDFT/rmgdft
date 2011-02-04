/************************** SVN Revision Information **************************
 **    $Id: latgen.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/latgen.c *****
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
 *   void latgen (int *ibrav, REAL *celldm, REAL *A0I, REAL *A1I, REAL *A2I, 
 *                REAL *OMEGAI, int *flag)
 *   sets up the crystallographic vectors a0, a1, and a2.
 * INPUTS
 *   ibrav: bravais lattice type
 *   celldm:  see the table below
 *   flag: if it is true, A0I, A1I, A2I are cell related (0-1.0)
 * OUTPUT
 *   A0I, A1I, A2I: lattice vector
 *   OMEGAI:  volume of a unit cell
 * PARENTS
 *   init.c
 * CHILDREN
 *   cross_product.c
 * SEE ALSO
 *   ibrav and celldm are defined in the table below
 *  Lattice constants  = celldm(1)
 * 
 * -----------------------------------------------------------------------
 * 
 * point group bravais lattice    ibrav  celldm(2)-celldm(6)
 * .......................................................................
 * 432,<4>3m,m3m     sc          1     not used
 * .......................................................................
 * 23,m3             sc          1         "
 * .......................................................................
 * 432,<4>3m,m3m    fcc          2         "
 * .......................................................................
 * 23,m3            fcc          2         "
 * .......................................................................
 * 432,<4>3m,m3m    bcc          3         "
 * .......................................................................
 * 23,m3            bcc          3         "
 * .......................................................................
 * 622,6mm,                            
 * <6>m2,6/mmm      hex(p)       4      celldm(3)=c/a
 * .......................................................................
 * 6,<6>,6/m,       hex(p)
 * 32,3m,<3>m      trig(p)       4         "
 * .......................................................................
 * 3,<3>           trig(p)       4         "
 * .......................................................................
 * 32,3m,<3>m      trig(r)       5     celldm(4)=cos(aalpha)
 * .......................................................................
 * 3,<3>           trig(r)       5         "
 * .......................................................................
 * 422,4mm, 
 * <4>2m,4/mmm     tetr(p)       6      celldm(3)=c/a
 * .......................................................................
 * 4,<4>,4/m       tetr(p)       6         "
 * .......................................................................
 * 422,4mm,
 * <4>2m,4/mmm     tetr(i)       7         "
 * .......................................................................
 * 4,<4>,4/m       tetr(i)       7         "
 * .......................................................................
 * 222,mm2,mmm     orth(p)       8     above + celldm(2)=b/a
 * .......................................................................
 * 2,m,2/m         mcln(p)      12     above + celldm(4)=cos(ab)
 * .......................................................................
 * 1,<1>           tcln(p)      14     celldm(2)= b/a
 *                                     celldm(3)= c/a
 *				       celldm(4)= cos(bc)
 *	  			       celldm(5)= cos(ac)
 *				       celldm(6)= cos(ab)
 * -----------------------------------------------------------------------
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"
#include "recips.h"

#define     SQRT2         1.414213562373
#define     SQRT3         1.732050807569


/* If flag is true then A0I,A1I,A2I are cell relative (0.0-1.0) */
void latgen (int *ibrav, REAL * celldm, REAL * A0I, REAL * A1I, REAL * A2I,
             REAL * OMEGAI, int *flag)
{

    int ir;
    REAL term, term1, term2, cbya, sine, singam;
    REAL distance, t1;
    REAL cvec[3];

    /* Initialise the appropriate variables */

    alat = celldm[0];

    for (ir = 0; ir < 3; ir++)
    {
        A0I[ir] = 0.0;
        A1I[ir] = 0.0;
        A2I[ir] = 0.0;
    }

    switch (*ibrav)
    {
    case CUBIC_PRIMITIVE:

        A0I[0] = celldm[0];
        A1I[1] = celldm[1];
        if (celldm[1] <= 0.0)
            A1I[1] = celldm[0];
        A2I[2] = celldm[2];
        if (celldm[2] <= 0.0)
            A2I[2] = celldm[0];
        break;

    case CUBIC_FC:

        term = alat / 2.0;
        A0I[0] = term;
        A0I[1] = term;
        A1I[1] = term;
        A1I[2] = term;
        A2I[0] = term;
        A2I[2] = term;
        break;

    case CUBIC_BC:

        term = alat / 2.0;
        for (ir = 0; ir < 3; ir++)
        {
            A0I[ir] = term;
            A1I[ir] = term;
            A2I[ir] = term;
        }                       /* for ir */
        A0I[2] = -term;
        A1I[0] = -term;
        A2I[1] = -term;
        break;

    case HEXAGONAL:
        cbya = celldm[2];
        A0I[0] = alat;
        A1I[0] = alat / 2.0;
        A1I[1] = alat * SQRT3 / 2.0;
        A2I[2] = alat * cbya;
        break;

    case TRIGONAL_PRIMITIVE:

        term1 = sqrt (1.0 + 2.0 * celldm[3]);
        term2 = sqrt (1.0 - celldm[3]);
        A0I[1] = SQRT2 * alat * term2 / SQRT3;
        A0I[2] = alat * term1 / SQRT3;
        A1I[0] = alat * term2 / SQRT2;
        A1I[1] = -A1I[0] / SQRT3;
        A1I[2] = A0I[2];
        A2I[0] = -A1I[0];
        A2I[1] = A1I[1];
        A2I[2] = A0I[2];
        break;

    case TETRAGONAL_PRIMITIVE:

        cbya = celldm[2];
        A0I[0] = alat;
        A1I[1] = alat;
        A2I[2] = alat * cbya;
        break;

    case TETRAGONAL_BC:

        cbya = celldm[2];
        A0I[0] = alat / 2.0;
        A0I[1] = A0I[0];
        A0I[2] = cbya * alat / 2.0;
        A1I[0] = A0I[0];
        A1I[1] = -A0I[0];
        A1I[2] = A0I[2];
        A2I[0] = -A0I[0];
        A2I[1] = -A0I[0];
        A2I[2] = A0I[2];
        break;

    case ORTHORHOMBIC_PRIMITIVE:

        A0I[0] = alat;
        A1I[1] = alat * celldm[1];
        A2I[2] = alat * celldm[2];
        break;

    case ORTHORHOMBIC_BASE_CENTRED:

        /* not programmed */
        error_handler ("bravais lattice not programmed.");
        break;

    case ORTHORHOMBIC_BC:

        /* not programmed */
        error_handler ("bravais lattice not programmed.");
        break;

    case ORTHORHOMBIC_FC:

        /* not programmed */
        error_handler ("bravais lattice not programmed.");
        break;

    case MONOCLINIC_PRIMITIVE:

        sine = sqrt (1.0 - celldm[3] * celldm[3]);
        A0I[0] = alat;
        A1I[0] = alat * celldm[1] * celldm[3];
        A1I[1] = alat * celldm[1] * sine;
        A2I[2] = alat * celldm[2];
        break;

    case MONOCLINIC_BASE_CENTRED:
        /* not programmed */
        error_handler ("bravais lattice not programmed.");
        break;

    case TRICLINIC_PRIMITIVE:

        singam = sqrt (1.0 - celldm[5] * celldm[5]);
        term = sqrt ((1.0 + 2.0 * celldm[3] * celldm[4] * celldm[5] -
                      celldm[3] * celldm[3]
                      - celldm[4] * celldm[4]
                      - celldm[5] * celldm[5]) / (1.0 - celldm[5] * celldm[5]));
        A0I[0] = alat;
        A1I[0] = alat * celldm[1] * celldm[5];
        A1I[1] = alat * celldm[1] * singam;
        A2I[0] = alat * celldm[2] * celldm[4];
        A2I[1] = alat * celldm[2] * (celldm[3] - celldm[4] * celldm[5]) / singam;
        A2I[2] = alat * celldm[2] * term;

        break;

    default:

        error_handler ("bravais lattice not programmed.");


    }                           /* end switch (*ibrav) */


    cross_product (A0I, A1I, cvec);
    *OMEGAI = cvec[0] * A2I[0] + cvec[1] * A2I[1] + cvec[2] * A2I[2];

    *OMEGAI = fabs (*OMEGAI);

    /* Generate volume element */
    t1 = (REAL) (ct.psi_nbasis);
    ct.vel = *OMEGAI / t1;

    t1 = (REAL) (ct.vh_nbasis);
    ct.vel_f = *OMEGAI / t1;

    /* Calculate length of supercell */
    distance = 0.0;
    for (ir = 0; ir < 3; ir++)
        distance += A0I[ir] * A0I[ir];

    ct.xside = sqrt (distance);

    distance = 0.0;
    for (ir = 0; ir < 3; ir++)
        distance += A1I[ir] * A1I[ir];

    ct.yside = sqrt (distance);

    distance = 0.0;
    for (ir = 0; ir < 3; ir++)
        distance += A2I[ir] * A2I[ir];

    ct.zside = sqrt (distance);

    /* Calculate grid size in crystal coordinates */

    t1 = (REAL) ct.psi_nxgrid;
    ct.hxgrid = 1.0 / t1;
    ct.hxxgrid = 1.0 / (REAL) ct.vh_nxgrid;

    t1 = (REAL) ct.psi_nygrid;
    ct.hygrid = 1.0 / t1;
    ct.hyygrid = 1.0 / (REAL) ct.vh_nygrid;

    t1 = (REAL) ct.psi_nzgrid;
    ct.hzgrid = 1.0 / t1;
    ct.hzzgrid = 1.0 / (REAL) ct.vh_nzgrid;

    if (*flag)
    {
        for (ir = 0; ir < 3; ir++)
        {
            A0I[ir] /= celldm[0];
            A1I[ir] /= celldm[0];
            A2I[ir] /= celldm[0];
        }                       /* end for */
    }                           /* end if */

}                               /* end latgen */

/******/
