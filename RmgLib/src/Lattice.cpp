/*
 *
 * Copyright (c) 1995, Emil Briggs
 * Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                     Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 * Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                     Marco Buongiorno Nardelli,Charles Brabec, 
 *                     Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                     Jerzy Bernholc
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * FUNCTION
 *   void latgen (int *ibrav, double *celldm, double *A0I, double *A1I, double *A2I, 
 *                double *OMEGAI, int *flag)
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
#include "Lattice.h"
#include "rmg_error.h"

#define         SQRT2       1.414213562373
#define         SQRT3       1.732050807569



void Lattice::cross_product (double * a, double * b, double * c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = -(a[0] * b[2] - a[2] * b[0]);
    c[2] = a[0] * b[1] - a[1] * b[0];
}                               /* end cross_product */

void Lattice::to_cartesian (double *crystal, double *cartesian)
{
    int ir;

    /*  position is in crystal coordinates, get cartesian coordinates */

    for (ir = 0; ir < 3; ir++)
    {
        cartesian[ir] = crystal[0] * Lattice::a0[ir] + crystal[1] * Lattice::a1[ir] + crystal[2] * Lattice::a2[ir];
    }                           /* end for ir */

}                               /* end to_cartesian */

void Lattice::to_crystal (double *crystal, double *cartesian)
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
void Lattice::latgen (double * celldm, double * OMEGAI, double *a0, double *a1, double *a2, int *flag)
{

    int ir;
    double term, term1, term2, cbya, sine, singam, alat;
    double distance;
    double cvec[3];

    /* Initialise the appropriate variables */

    alat = celldm[0];

    for (ir = 0; ir < 3; ir++)
    {
        Lattice::a0[ir] = 0.0;
        Lattice::a1[ir] = 0.0;
        Lattice::a2[ir] = 0.0;
        Lattice::celldm[ir] = celldm[ir];
    }

    switch (Lattice::ibrav)
    {


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

        case CUBIC_PRIMITIVE:
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
    Lattice::omega = *OMEGAI;

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

double Lattice::metric (double * crystal)
{
    double cartesian[3];          /* cartesian coordinates of point */
    double distance;
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
double Lattice::get_omega(void)
{
    return Lattice::omega;
}
double Lattice::get_celldm(int which)
{
    return Lattice::celldm[which];
}

double Lattice::get_a0(int which)
{
    return Lattice::a0[which];
}
double Lattice::get_a1(int which)
{
    return Lattice::a1[which];
}
double Lattice::get_a2(int which)
{
    return Lattice::a2[which];
}
double Lattice::get_b0(int which)
{
    return Lattice::b0[which];
}
double Lattice::get_b1(int which)
{
    return Lattice::b1[which];
}
double Lattice::get_b2(int which)
{
    return Lattice::b2[which];
}

double Lattice::get_xside(void)
{
    return Lattice::xside;
}
double Lattice::get_yside(void)
{
    return Lattice::yside;
}
double Lattice::get_zside(void)
{
    return Lattice::zside;
}

//****/

