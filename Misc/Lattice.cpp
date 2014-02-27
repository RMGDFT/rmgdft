/************************** SVN Revision Information **************************
 **    $Id: latgen.c 2099 2014-02-09 16:44:35Z ebriggs $    **
******************************************************************************/

/**f* QMD-MGDFT/latgen.c *****
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
 *   void latgen (int *ibrav, rmg_double_t *celldm, rmg_double_t *A0I, rmg_double_t *A1I, rmg_double_t *A2I, 
 *                rmg_double_t *OMEGAI, int *flag)
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

#include <math.h>
#include "BaseGrid.h"
#include "Lattice.h"
#include "rmg_error.h"



using namespace std;


void Lattice::cross_product (rmg_double_t * a, rmg_double_t * b, rmg_double_t * c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = -(a[0] * b[2] - a[2] * b[0]);
    c[2] = a[0] * b[1] - a[1] * b[0];
}                               /* end cross_product */

void Lattice::to_cartesian (rmg_double_t *crystal, rmg_double_t *cartesian)
{
    int ir;

    /*  position is in crystal coordinates, get cartesian coordinates */

    for (ir = 0; ir < 3; ir++)
    {
        cartesian[ir] = crystal[0] * Lattice::a0[ir] + crystal[1] * Lattice::a1[ir] + crystal[2] * Lattice::a2[ir];
    }                           /* end for ir */

}                               /* end to_cartesian */

void Lattice::to_crystal (rmg_double_t *crystal, rmg_double_t *cartesian)
{
    int ir;
    BaseGrid G;

    if (G.ibrav == HEXAGONAL)
    {

        crystal[0] = (cartesian[0] - cartesian[1] / SQRT3) / Lattice::celldm[0];
        crystal[1] = cartesian[1] / (SQRT3 / 2.0) / Lattice::celldm[0];
        crystal[2] = cartesian[2] * b2[2];

        if (crystal[0] < 0.0)
            crystal[0] += 1.0;
        if (crystal[1] < 0.0)
            crystal[1] += 1.0;
        if (crystal[2] < 0.0)
            crystal[2] += 1.0;
        if (crystal[0] > 1.0)
            crystal[0] -= 1.0;
        if (crystal[1] > 1.0)
            crystal[1] -= 1.0;
        if (crystal[2] > 1.0)
            crystal[2] -= 1.0;

    }
    else if (G.ibrav == CUBIC_PRIMITIVE)
    {

        crystal[0] = cartesian[0] / Lattice::celldm[0];
        crystal[1] = cartesian[1] / Lattice::celldm[0];
        crystal[2] = cartesian[2] / Lattice::celldm[0];

        if (crystal[0] < 0.0)
            crystal[0] += 1.0;
        if (crystal[1] < 0.0)
            crystal[1] += 1.0;
        if (crystal[2] < 0.0)
            crystal[2] += 1.0;
        if (crystal[0] > 1.0)
            crystal[0] -= 1.0;
        if (crystal[1] > 1.0)
            crystal[1] -= 1.0;
        if (crystal[2] > 1.0)
            crystal[2] -= 1.0;

    }
    else if (G.ibrav == ORTHORHOMBIC_PRIMITIVE)
    {

        crystal[0] = cartesian[0] / Lattice::celldm[0];
        crystal[1] = cartesian[1] / (Lattice::celldm[0] * Lattice::celldm[1]);
        crystal[2] = cartesian[2] / (Lattice::celldm[0] * Lattice::celldm[2]);

        if (crystal[0] < 0.0)
            crystal[0] += 1.0;
        if (crystal[1] < 0.0)
            crystal[1] += 1.0;
        if (crystal[2] < 0.0)
            crystal[2] += 1.0;
        if (crystal[0] > 1.0)
            crystal[0] -= 1.0;
        if (crystal[1] > 1.0)
            crystal[1] -= 1.0;
        if (crystal[2] > 1.0)
            crystal[2] -= 1.0;

    }
    else
    {

        for (ir = 0; ir < 3; ir++)
        {
            crystal[ir] = cartesian[0] * b0[ir] + cartesian[1] * b1[ir] + cartesian[2] * b2[ir];
        }                       /* end for ir */

    }                           /* end if */

}                               /* end to_crystal  */

void Lattice::recips (void)
{

    double mag;

    cross_product (Lattice::a1, Lattice::a2, Lattice::b0);
    mag = Lattice::a0[0] * Lattice::b0[0] + Lattice::a0[1] * Lattice::b0[1] + Lattice::a0[2] * Lattice::b0[2];
    Lattice::b0[0] /= mag;
    Lattice::b0[1] /= mag;
    Lattice::b0[2] /= mag;
    cross_product (Lattice::a2, Lattice::a0, Lattice::b1);
    Lattice::b1[0] /= mag;
    Lattice::b1[1] /= mag;
    Lattice::b1[2] /= mag;
    cross_product (Lattice::a0, Lattice::a1, Lattice::b2);
    Lattice::b2[0] /= mag;
    Lattice::b2[1] /= mag;
    Lattice::b2[2] /= mag;


}                               /* end recips */


// If flag is true then A0I,A1I,A2I are cell relative (0.0-1.0)
void Lattice::latgen (rmg_double_t * celldm, rmg_double_t * OMEGAI, int *flag)
{

    int ir;
    rmg_double_t term, term1, term2, cbya, sine, singam, alat;
    rmg_double_t distance, t1;
    rmg_double_t cvec[3];
    BaseGrid G;

    /* Initialise the appropriate variables */

    alat = celldm[0];

    for (ir = 0; ir < 3; ir++)
    {
        Lattice::a0[ir] = 0.0;
        Lattice::a1[ir] = 0.0;
        Lattice::a2[ir] = 0.0;
    }

    switch (G.ibrav)
    {

        case CUBIC_PRIMITIVE:
            Lattice::a0[0] = celldm[0];
            Lattice::a1[1] = celldm[1];
            if (celldm[1] <= 0.0)
                Lattice::a1[1] = celldm[0];
            Lattice::a2[2] = celldm[2];
            if (celldm[2] <= 0.0)
                Lattice::a2[2] = celldm[0];
            break;

        case CUBIC_FC:

            term = alat / 2.0;
            Lattice::a0[0] = term;
            Lattice::a0[1] = term;
            Lattice::a1[1] = term;
            Lattice::a1[2] = term;
            Lattice::a2[0] = term;
            Lattice::a2[2] = term;
            break;

        case CUBIC_BC:

            term = alat / 2.0;
            for (ir = 0; ir < 3; ir++)
            {
                Lattice::a0[ir] = term;
                Lattice::a1[ir] = term;
                Lattice::a2[ir] = term;
            }                       /* for ir */
            Lattice::a0[2] = -term;
            Lattice::a1[0] = -term;
            Lattice::a2[1] = -term;
            break;

        case HEXAGONAL:
            cbya = celldm[2];
            Lattice::a0[0] = alat;
            Lattice::a1[0] = alat / 2.0;
            Lattice::a1[1] = alat * SQRT3 / 2.0;
            Lattice::a2[2] = alat * cbya;
            break;

        case TRIGONAL_PRIMITIVE:

            term1 = sqrt (1.0 + 2.0 * celldm[3]);
            term2 = sqrt (1.0 - celldm[3]);
            Lattice::a0[1] = SQRT2 * alat * term2 / SQRT3;
            Lattice::a0[2] = alat * term1 / SQRT3;
            Lattice::a1[0] = alat * term2 / SQRT2;
            Lattice::a1[1] = -Lattice::a1[0] / SQRT3;
            Lattice::a1[2] = Lattice::a0[2];
            Lattice::a2[0] = -Lattice::a1[0];
            Lattice::a2[1] = Lattice::a1[1];
            Lattice::a2[2] = Lattice::a0[2];
            break;

        case TETRAGONAL_PRIMITIVE:

            cbya = celldm[2];
            Lattice::a0[0] = alat;
            Lattice::a1[1] = alat;
            Lattice::a2[2] = alat * cbya;
            break;

        case TETRAGONAL_BC:

            cbya = celldm[2];
            Lattice::a0[0] = alat / 2.0;
            Lattice::a0[1] = Lattice::a0[0];
            Lattice::a0[2] = cbya * alat / 2.0;
            Lattice::a1[0] = Lattice::a0[0];
            Lattice::a1[1] = -Lattice::a0[0];
            Lattice::a1[2] = Lattice::a0[2];
            Lattice::a2[0] = -Lattice::a0[0];
            Lattice::a2[1] = -Lattice::a0[0];
            Lattice::a2[2] = Lattice::a0[2];
            break;

        case ORTHORHOMBIC_PRIMITIVE:

            Lattice::a0[0] = alat;
            Lattice::a1[1] = alat * celldm[1];
            Lattice::a2[2] = alat * celldm[2];
            break;

        case ORTHORHOMBIC_BASE_CENTRED:

            /* not programmed */
            rmg_error_handler ("bravais lattice not programmed.");
            break;

        case ORTHORHOMBIC_BC:

            /* not programmed */
            rmg_error_handler ("bravais lattice not programmed.");
            break;

        case ORTHORHOMBIC_FC:

            /* not programmed */
            rmg_error_handler ("bravais lattice not programmed.");
            break;

        case MONOCLINIC_PRIMITIVE:

            sine = sqrt (1.0 - celldm[3] * celldm[3]);
            Lattice::a0[0] = alat;
            Lattice::a1[0] = alat * celldm[1] * celldm[3];
            Lattice::a1[1] = alat * celldm[1] * sine;
            Lattice::a2[2] = alat * celldm[2];
            break;

        case MONOCLINIC_BASE_CENTRED:
            /* not programmed */
            rmg_error_handler ("bravais lattice not programmed.");
            break;

        case TRICLINIC_PRIMITIVE:

            singam = sqrt (1.0 - celldm[5] * celldm[5]);
            term = sqrt ((1.0 + 2.0 * celldm[3] * celldm[4] * celldm[5] -
                          celldm[3] * celldm[3]
                          - celldm[4] * celldm[4]
                          - celldm[5] * celldm[5]) / (1.0 - celldm[5] * celldm[5]));
            Lattice::a0[0] = alat;
            Lattice::a1[0] = alat * celldm[1] * celldm[5];
            Lattice::a1[1] = alat * celldm[1] * singam;
            Lattice::a2[0] = alat * celldm[2] * celldm[4];
            Lattice::a2[1] = alat * celldm[2] * (celldm[3] - celldm[4] * celldm[5]) / singam;
            Lattice::a2[2] = alat * celldm[2] * term;

            break;

        default:

            printf ("ibrav is set to %d", G.ibrav);
            rmg_error_handler ("bravais lattice not programmed.");


    }                           /* end switch (*ibrav) */

    cross_product (Lattice::a0, Lattice::a1, cvec);
    *OMEGAI = cvec[0] * Lattice::a2[0] + cvec[1] * Lattice::a2[1] + cvec[2] * Lattice::a2[2];

    *OMEGAI = fabs (*OMEGAI);

    /* Generate volume element */
    t1 = (rmg_double_t) (G.NX_GRID * G.NY_GRID * G.NZ_GRID);
    Lattice::vel = *OMEGAI / t1;

    t1 = (rmg_double_t) (G.FNX_GRID * G.FNY_GRID * G.FNZ_GRID);
    Lattice::vel_f = *OMEGAI / t1;

    /* Calculate length of supercell */
    distance = 0.0;
    for (ir = 0; ir < 3; ir++)
        distance += Lattice::a0[ir] * Lattice::a0[ir];

    Lattice::xside = sqrt (distance);

    distance = 0.0;
    for (ir = 0; ir < 3; ir++)
        distance += Lattice::a1[ir] * Lattice::a1[ir];

    Lattice::yside = sqrt (distance);

    distance = 0.0;
    for (ir = 0; ir < 3; ir++)
        distance += Lattice::a2[ir] * Lattice::a2[ir];

    Lattice::zside = sqrt (distance);

    /* Calculate grid size in crystal coordinates */

    t1 = (rmg_double_t) G.NX_GRID;
    Lattice::hxgrid = 1.0 / t1;
    Lattice::hxxgrid = Lattice::hxgrid / (rmg_double_t) G.FG_NX;

    t1 = (rmg_double_t) G.NY_GRID;
    Lattice::hygrid = 1.0 / t1;
    Lattice::hyygrid = Lattice::hygrid / (rmg_double_t) G.FG_NY;

    t1 = (rmg_double_t) G.NZ_GRID;
    Lattice::hzgrid = 1.0 / t1;
    Lattice::hzzgrid = Lattice::hzgrid / (rmg_double_t) G.FG_NZ;

    if (*flag)
    {
        for (ir = 0; ir < 3; ir++)
        {
            Lattice::a0[ir] /= celldm[0];
            Lattice::a1[ir] /= celldm[0];
            Lattice::a2[ir] /= celldm[0];
        }                       /* end for */
    }                           /* end if */

    Lattice::recips();
}                               /* end latgen */

rmg_double_t Lattice::metric (rmg_double_t * crystal)
{
    rmg_double_t cartesian[3];          /* cartesian coordinates of point */
    rmg_double_t distance;
    int ir;
    Lattice::to_cartesian (crystal, cartesian);

    distance = 0.0;

    for (ir = 0; ir < 3; ir++)
        distance += cartesian[ir] * cartesian[ir];

    distance = sqrt (distance);

    return (distance);

}                               /* end metric */

//****/

// lengths of the sides of the supercell
rmg_double_t Lattice::xside;
rmg_double_t Lattice::yside;
rmg_double_t Lattice::zside;

// lattice vectors
rmg_double_t Lattice::a0[3];
rmg_double_t Lattice::a1[3];
rmg_double_t Lattice::a2[3];

// reciprocal lattice vectors
rmg_double_t Lattice::b0[3];
rmg_double_t Lattice::b1[3];
rmg_double_t Lattice::b2[3];      

// cell dimensions
rmg_double_t Lattice::celldm[6];

// Total cell volume */
rmg_double_t Lattice::omega;

// Volume elements on coarse and fine grids
rmg_double_t Lattice::vel;
rmg_double_t Lattice::vel_f;

// Global uniform grid spacing in x
rmg_double_t Lattice::hxgrid;

// Global uniform grid spacing in y
rmg_double_t Lattice::hygrid;

// Global uniform grid spacing in z
rmg_double_t Lattice::hzgrid;

// The fine uniform grid spacing in x
rmg_double_t Lattice::hxxgrid;

// The fine uniform grid spacing in y
rmg_double_t Lattice::hyygrid;

// The fine uniform grid spacing in z
rmg_double_t Lattice::hzzgrid;


extern "C" rmg_double_t get_xside(void)
{
    return Lattice::xside;
}
extern "C" rmg_double_t get_yside(void)
{
    return Lattice::yside;
}
extern "C" rmg_double_t get_zside(void)
{
    return Lattice::zside;
}
extern "C" rmg_double_t get_hxgrid(void)
{
    return Lattice::hxgrid;
}
extern "C" rmg_double_t get_hygrid(void)
{
    return Lattice::hygrid;
}
extern "C" rmg_double_t get_hzgrid(void)
{
    return Lattice::hzgrid;
}
extern "C" rmg_double_t get_hxxgrid(void)
{
    return Lattice::hxxgrid;
}
extern "C" rmg_double_t get_hyygrid(void)
{
    return Lattice::hyygrid;
}
extern "C" rmg_double_t get_hzzgrid(void)
{
    return Lattice::hzzgrid;
}
extern "C" rmg_double_t get_vel(void)
{
    return Lattice::vel;
}
extern "C" rmg_double_t get_vel_f(void)
{
    return Lattice::vel_f;
}
extern "C" rmg_double_t get_celldm(int which)
{
    return Lattice::celldm[which];
}
extern "C" rmg_double_t get_a0(int which)
{
    return Lattice::a0[which];
}
extern "C" rmg_double_t get_a1(int which)
{
    return Lattice::a1[which];
}
extern "C" rmg_double_t get_a2(int which)
{
    return Lattice::a2[which];
}
extern "C" rmg_double_t get_b0(int which)
{
    return Lattice::b0[which];
}
extern "C" rmg_double_t get_b1(int which)
{
    return Lattice::b1[which];
}
extern "C" rmg_double_t get_b2(int which)
{
    return Lattice::b2[which];
}
extern "C" void to_crystal (rmg_double_t *crystal, rmg_double_t *cartesian)
{
    Lattice L;
    L.to_crystal(crystal, cartesian);
}
extern "C" void to_cartesian (rmg_double_t *crystal, rmg_double_t *cartesian)
{
    Lattice L;
    L.to_cartesian(crystal, cartesian);
}
extern "C" void cross_product (rmg_double_t * a, rmg_double_t * b, rmg_double_t * c)
{
    Lattice L;
    L.cross_product(a, b, c);
}
extern "C" rmg_double_t metric (rmg_double_t * crystal)
{
    Lattice L;
    return L.metric(crystal);
}
extern "C" void recips(void)
{
    Lattice L;
    L.recips();
}
extern "C" void latgen (rmg_double_t *celldm, rmg_double_t *a0, rmg_double_t *a1, rmg_double_t *a2, rmg_double_t *OMEGAI, int *flag)
{
    Lattice L;
    L.latgen(celldm, OMEGAI, flag);
    a0[0] = Lattice::a0[0];
    a0[1] = Lattice::a0[1];
    a0[2] = Lattice::a0[2];
    a1[0] = Lattice::a1[0];
    a1[1] = Lattice::a1[1];
    a1[2] = Lattice::a1[2];
    a2[0] = Lattice::a2[0];
    a2[1] = Lattice::a2[1];
    a2[2] = Lattice::a2[2];

}

