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

#define         SQRT2       1.414213562373
#define         SQRT3       1.732050807569


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

    if (Lattice::ibrav == HEXAGONAL)
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
    else if (Lattice::ibrav == CUBIC_PRIMITIVE)
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
    else if (Lattice::ibrav == ORTHORHOMBIC_PRIMITIVE)
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
void Lattice::latgen (rmg_double_t * celldm, rmg_double_t * OMEGAI, rmg_double_t *a0, rmg_double_t *a1, rmg_double_t *a2, int *flag)
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

    switch (Lattice::ibrav)
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
            rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
            break;

        case ORTHORHOMBIC_BC:

            /* not programmed */
            rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
            break;

        case ORTHORHOMBIC_FC:

            /* not programmed */
            rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
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
            rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
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

            printf ("ibrav is set to %d", Lattice::ibrav);
            rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");


    }                           /* end switch (*ibrav) */

    cross_product (Lattice::a0, Lattice::a1, cvec);
    *OMEGAI = cvec[0] * Lattice::a2[0] + cvec[1] * Lattice::a2[1] + cvec[2] * Lattice::a2[2];

    *OMEGAI = fabs (*OMEGAI);

    /* Generate volume element */
    t1 = (rmg_double_t) (G.get_NX_GRID() * G.get_NY_GRID() * G.get_NZ_GRID());
    Lattice::vel = *OMEGAI / t1;

    t1 = (rmg_double_t) (G.get_FNX_GRID() * G.get_FNY_GRID() * G.get_FNZ_GRID());
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

    a0[0] = Lattice::a0[0];
    a0[1] = Lattice::a0[1];
    a0[2] = Lattice::a0[2];
    a1[0] = Lattice::a1[0];
    a1[1] = Lattice::a1[1];
    a1[2] = Lattice::a1[2];
    a2[0] = Lattice::a2[0];
    a2[1] = Lattice::a2[1];
    a2[2] = Lattice::a2[2];

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


int Lattice::get_ibrav_type(void)
{
    return Lattice::ibrav;
}
void Lattice::set_ibrav_type(int newtype)
{
  Lattice::ibrav = newtype;
}
rmg_double_t Lattice::get_xside(void)
{
    return Lattice::xside;
}
rmg_double_t Lattice::get_yside(void)
{
    return Lattice::yside;
}
rmg_double_t Lattice::get_zside(void)
{
    return Lattice::zside;
}

//****/

// Grid bravais lattice type 
int Lattice::ibrav;

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

/// C interface function
extern "C" void set_ibrav_type(int newtype)
{
  Lattice L;
  L.set_ibrav_type(newtype);
}
/// C interface function
extern "C" int get_ibrav_type(void)
{
  Lattice L;
  return L.get_ibrav_type();
}
extern "C" rmg_double_t get_xside(void)
{
    Lattice L;
    return L.get_xside();
}
extern "C" rmg_double_t get_yside(void)
{
    Lattice L;
    return L.get_yside();
}
extern "C" rmg_double_t get_zside(void)
{
    Lattice L;
    return L.get_zside();
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
    L.latgen(celldm, OMEGAI, a0, a1, a2, flag);
}

