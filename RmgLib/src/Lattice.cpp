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

#include <cmath>
#include <boost/multi_array.hpp>
#include "Lattice.h"
#include "rmg_error.h"
#include "blas.h"

//#define         SQRT2       1.414213562373
//#define         SQRT3       1.732050807569


double Lattice::dot_product(double *a, double *b)
{
    double ret = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return ret;
}
void Lattice::cross_product (double * a, double * b, double * c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = -(a[0] * b[2] - a[2] * b[0]);
    c[2] = a[0] * b[1] - a[1] * b[0];
}                               /* end cross_product */

void Lattice::to_cartesian (double *crystal, double *cartesian)
{

    /*  position is in crystal coordinates, get cartesian coordinates */

    for (int ir = 0; ir < 3; ir++)
    {
        cartesian[ir] = crystal[0] * Lattice::a0[ir] + crystal[1] * Lattice::a1[ir] + crystal[2] * Lattice::a2[ir];
    }                           /* end for ir */

}                               /* end to_cartesian */

void Lattice::to_cartesian_input (double *crystal, double *cartesian)
{

    /*  position is in crystal coordinates, get cartesian coordinates using the input (non-rotated lattice vectors) */

    for (int ir = 0; ir < 3; ir++)
    {
        cartesian[ir] = crystal[0] * Lattice::a0i[ir] + crystal[1] * Lattice::a1i[ir] + crystal[2] * Lattice::a2i[ir];
    }                           /* end for ir */

}                               /* end to_cartesian */

void Lattice::to_crystal (double *crystal, double *cartesian)
{

    crystal[0] = cartesian[0] * b0[0] + cartesian[1] * b0[1] + cartesian[2] * b0[2];
    crystal[1] = cartesian[0] * b1[0] + cartesian[1] * b1[1] + cartesian[2] * b1[2];
    crystal[2] = cartesian[0] * b2[0] + cartesian[1] * b2[1] + cartesian[2] * b2[2];

    if (crystal[0] < -1.0e-10)
        crystal[0] += 1.0;
    if (crystal[1] < -1.0e-10)
        crystal[1] += 1.0;
    if (crystal[2] < -1.0e-10)
        crystal[2] += 1.0;
    if (crystal[0] > 1.0)
        crystal[0] -= 1.0;
    if (crystal[1] > 1.0)
        crystal[1] -= 1.0;
    if (crystal[2] > 1.0)
        crystal[2] -= 1.0;


}                               /* end to_crystal  */

// Transforms vector
void Lattice::to_crystal_vector (double *crystal, double *cartesian)
{

    crystal[0] = cartesian[0] * b0[0] + cartesian[1] * b0[1] + cartesian[2] * b0[2];
    crystal[1] = cartesian[0] * b1[0] + cartesian[1] * b1[1] + cartesian[2] * b1[2];
    crystal[2] = cartesian[0] * b2[0] + cartesian[1] * b2[1] + cartesian[2] * b2[2];
}

void Lattice::to_crystal_half (double *crystal, double *cartesian)
{

    crystal[0] = cartesian[0] * b0[0] + cartesian[1] * b0[1] + cartesian[2] * b0[2];
    crystal[1] = cartesian[0] * b1[0] + cartesian[1] * b1[1] + cartesian[2] * b1[2];
    crystal[2] = cartesian[0] * b2[0] + cartesian[1] * b2[1] + cartesian[2] * b2[2];

    if (crystal[0] < -0.5)
        crystal[0] += 1.0;
    if (crystal[1] < -0.5)
        crystal[1] += 1.0;
    if (crystal[2] < -0.5)
        crystal[2] += 1.0;
    if (crystal[0] > 0.5)
        crystal[0] -= 1.0;
    if (crystal[1] > 0.5)
        crystal[1] -= 1.0;
    if (crystal[2] > 0.5)
        crystal[2] -= 1.0;


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


// flag == true: determine celldm from a0, a1, a2
// flag == false: determin a0, a1, a2 from celldm
void Lattice::latgen (double * celldm, double * OMEGAI, double *a0, double *a1, double *a2, bool flag)
{

    int ir;
    double term, term1, term2, cbya, sine, singam, alat;
    double distance;
    double cvec[3];

    /* Initialise the appropriate variables */

    if(Lattice::ibrav == None || flag)
    {
        for (ir = 0; ir < 3; ir++)
        {
            Lattice::a0[ir] = a0[ir];
            Lattice::a1[ir] = a1[ir];
            Lattice::a2[ir] = a2[ir];
        }

        distance = 0.0;
        for (ir = 0; ir < 3; ir++)
            distance += Lattice::a0[ir] * Lattice::a0[ir];
        Lattice::celldm[0] = sqrt(distance);

        distance = 0.0;
        for (ir = 0; ir < 3; ir++)
            distance += Lattice::a1[ir] * Lattice::a1[ir];
        Lattice::celldm[1] = sqrt(distance)/Lattice::celldm[0];

        distance = 0.0;
        for (ir = 0; ir < 3; ir++)
            distance += Lattice::a2[ir] * Lattice::a2[ir];
        Lattice::celldm[2] = sqrt(distance)/Lattice::celldm[0];


        alat = Lattice::celldm[0];

        for (ir = 0; ir < 6; ir++)
            celldm[ir] = Lattice::celldm[ir];

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

        Lattice::recips();
        return;
    }
    else
    {
        for (ir = 0; ir < 6; ir++)
            Lattice::celldm[ir] = celldm[ir];
        // force the lattice symmetry, lattice vectors are determined from celldm[0],[1],[2]
        for (ir = 0; ir < 3; ir++)
        {
            Lattice::a0[ir] = 0.0;
            Lattice::a1[ir] = 0.0;
            Lattice::a2[ir] = 0.0;
        }
    }
    if(Lattice::ibrav == None) return;

    alat = celldm[0];
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
            Lattice::a1[0] = -alat / 2.0;
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

            if(celldm[1] <= 0.0) rmg_error_handler (__FILE__, __LINE__, "celldm[1] must be positive.");
            if(celldm[2] <= 0.0) rmg_error_handler (__FILE__, __LINE__, "celldm[2] must be positive.");
            
            Lattice::a0[0] = 0.5 * celldm[0];
            Lattice::a0[1] = Lattice::a0[0] * celldm[1];
            Lattice::a0[2] = Lattice::a0[0] * celldm[2];
            Lattice::a1[0] = - Lattice::a0[0];
            Lattice::a1[1] = Lattice::a0[1];
            Lattice::a1[2] = Lattice::a0[2];
            Lattice::a2[0] = - Lattice::a0[0];
            Lattice::a2[1] = - Lattice::a0[1];
            Lattice::a2[2] = Lattice::a0[2];
            break;

        case ORTHORHOMBIC_FC:

            /* not programmed */
            rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
            break;

        case -MONOCLINIC_PRIMITIVE:
            Lattice::a0[0] = celldm[0];
            Lattice::a1[1] = celldm[0]*celldm[1];
            Lattice::a2[0] = celldm[0]*celldm[2]*celldm[4];
            Lattice::a2[2] = celldm[0]*celldm[2] * sqrt(1.0 - celldm[4]*celldm[4]);
            break;
 
        case MONOCLINIC_PRIMITIVE:

            //rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
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

            //rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
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

    Lattice::recips();

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

void Lattice::move_cell(double dt, int *cell_movable)
{
    for (int i = 0; i < 3; i++)
    {
        Lattice::a0[i] += dt * Lattice::cell_velocity[0*3 + i] * cell_movable[0*3+i];
        Lattice::a1[i] += dt * Lattice::cell_velocity[1*3 + i] * cell_movable[1*3+i];
        Lattice::a2[i] += dt * Lattice::cell_velocity[2*3 + i] * cell_movable[2*3+i];
    }

    double celldm[6]= {1.0,1.0,1.0,0.0,0.0,0.0},omega;
    Lattice::latgen (celldm, &omega, Lattice::a0, Lattice::a1, Lattice::a2, true);

}

bool zcheck(double x)
{
  return (std::abs(x) < 1.0e-4);
}
// Generates crystallographic parameters a,b,c from input lattice vectors
void Lattice::lat2abc(double *a0, double *a1, double *a2)
{
    a = sqrt(a0[0]*a0[0] + a0[1]*a0[1] + a0[2]*a0[2]);
    b = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    c = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);

    if(fabs(a-b) < 1.0e-4) b = a;
    if(fabs(a-c) < 1.0e-4) c = a;
    cosab = (a0[0]*a1[0] + a0[1]*a1[1] + a0[2]*a1[2])/a/b;
    cosac = (a0[0]*a2[0] + a0[1]*a2[1] + a0[2]*a2[2])/a/c;
    cosbc = (a2[0]*a1[0] + a2[1]*a1[1] + a2[2]*a1[2])/c/b;
    if(fabs(cosab) < 1.0e-4) cosab = 0.0;
    if(fabs(cosac) < 1.0e-4) cosac = 0.0;
    if(fabs(cosbc) < 1.0e-4) cosbc = 0.0;
    if(fabs(cosab - cosac) < 1.0e-4) cosac = cosab;
    if(fabs(cosab - cosbc) < 1.0e-4) cosbc = cosab;


    // Now figure out the lattice type
    if ( (zcheck(b) && zcheck(c)) || (zcheck(b-a) && zcheck(c-a)) )
    {
        // Case: a=b=c
        if ( zcheck(cosab) && zcheck(cosac) && zcheck(cosbc)) {
            // Case: simple cubic
            ibrav = 1;
        }
        else if ( !zcheck(cosab) && zcheck(cosac) && zcheck(cosbc)) {
            //Case: trigonal R with threefold axis along <111>
            ibrav =-5;
        }
        else if ( zcheck(cosab-0.5) && cosac == cosab && cosbc == cosab ) {
            //Case: fcc
            ibrav = 2;
        }
        else if ( std::abs ( cosab + 1.0/sqrt(3.0) ) < 1.0e-6 && cosac == cosab && cosbc == cosab ) {
            //Case: bcc with symmetric basis
            ibrav =-3;
        }
        else {
           //Case: unknown (possibly wrong)
           ibrav =-1;
        }
    }   
    else if ( (zcheck(b) && c > 0.0) || (b == a && c != a) )
    {
        //Case: a=b!=c
        if ( zcheck(cosab) && zcheck(cosac) && zcheck(cosbc)) {
            //Case: simple tetragonal
            ibrav = 6;
         }
         else if ( zcheck(cosab+0.5) && zcheck(cosac) && zcheck(cosbc) ) {
             //Case: simple hexagonal with 120 degree ab
             ibrav = 4;
         }
         else if ( zcheck(cosab-0.5) && zcheck(cosac) && zcheck(cosbc) ) {
             //Case: simple hexagonal with 60 degree ab
             ibrav = 15;
         }
         else
         {
             //Case: unknown (body-centered tetragonal or wrong data)
             ibrav =-1;
         }
    }   
    else if ( b > 0.0 && c > 0.0 && b != a && c != a )
    {
        //Case: a!=b!=c
        if ( zcheck(cosab) && zcheck(cosac) && zcheck(cosbc) ) {
            //Case: simple orthorhombic
            ibrav = 8;
        }
        else if ( !zcheck(cosab) && zcheck(cosac) && zcheck(cosbc) ) {
            //Case: monoclinic P, unique axis c
            ibrav = 12;
        }
        else if ( zcheck(cosab) && !zcheck(cosac) && zcheck(cosbc) ) {
            //Case: monoclinic P, unique axis b
            ibrav =-12;
        }
        else if ( !zcheck(cosab) && !zcheck(cosac) && !zcheck(cosbc) ) {
            //Case: triclinic
            ibrav = 14;
        }
        else { 
            //Case: unknown (base-, face-, body-centered orthorombic,
            //              (base-centered monoclinic)
            ibrav = -1;
        }
    }

//    if(ibrav < 0)
//        rmg_error_handler (__FILE__, __LINE__, "Negative ibrav not supported.\n");
}

void Lattice::abc2celldm(void)
{
#if 0
  IF (a <= 0.0_dp) CALL errore('abc2celldm','incorrect lattice parameter (a)',1)
  IF (b <  0.0_dp) CALL errore('abc2celldm','incorrect lattice parameter (b)',1)
  IF (c <  0.0_dp) CALL errore('abc2celldm','incorrect lattice parameter (c)',1)
  IF ( ABS (cosab) > 1.0_dp) CALL errore('abc2celldm', &
                   'incorrect lattice parameter (cosab)',1)
  IF ( ABS (cosac) > 1.0_dp) CALL errore('abc2celldm', &
                   'incorrect lattice parameter (cosac)',1)
  IF ( ABS (cosbc) > 1.0_dp) CALL errore('abc2celldm', &
       'incorrect lattice parameter (cosbc)',1)
#endif
  celldm[0] = a;
  celldm[1] = b / a;
  celldm[2] = c / a;

  if ( ibrav == 14 )
  {
     // ... triclinic lattice
     celldm[3] = cosbc;
     celldm[4] = cosac;
     celldm[5] = cosab;
  }
  else if ( ibrav ==-12 )
  {
     // ... monoclinic P lattice, unique axis b
     celldm[4] = cosac;
  }
  else
  {
     // ... trigonal and monoclinic lattices, unique axis c
     celldm[3] = cosab;
  }

}


static bool eqq(double x, double y)
{
    return (std::abs(x-y) < 1.0e-5);
}
static bool neqq(double x, double y)
{
    return (std::abs(x-y) >= 1.0e-5);
}

int Lattice::lat2ibrav (double *a0, double *a1, double *a2)
{
    //
    //     Returns ibrav from lattice vectors if recognized, 0 otherwise
    //
    a = sqrt(a0[0]*a0[0] + a0[1]*a0[1] + a0[2]*a0[2]);
    b = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    c = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);

    if(fabs(a-b) < 1.0e-5) b = a;
    if(fabs(a-c) < 1.0e-5) c = a;

    cosab = (a0[0]*a1[0] + a0[1]*a1[1] + a0[2]*a1[2])/a/b;
    cosac = (a0[0]*a2[0] + a0[1]*a2[1] + a0[2]*a2[2])/a/c;
    cosbc = (a2[0]*a1[0] + a2[1]*a1[1] + a2[2]*a1[2])/c/b;

    if(fabs(cosab) < 1.0e-5) cosab = 0.0;
    if(fabs(cosac) < 1.0e-5) cosac = 0.0;
    if(fabs(cosbc) < 1.0e-5) cosbc = 0.0;
    if(fabs(cosab - cosac) < 1.0e-5) cosac = cosab;
    if(fabs(cosab - cosbc) < 1.0e-5) cosbc = cosab;
    //
    // Assume triclinic if nothing suitable found
    //
    ibrav = 14;
    if ( eqq(a,b) && eqq(a,c) )
    {
        // Case: a=b=c
        if (eqq(cosab,cosac) && eqq(cosab,cosbc))
        {
            // Case: alpha = beta = gamma
            if ( eqq(cosab,0.0) )
            {
                // Cubic P - ibrav=1
                ibrav = 1;
            }
            else if ( eqq(cosab,0.5) )
            {
                // Cubic F - ibrav=2
                ibrav = 2;
            }
            else if ( eqq(cosab,-1.0/3.0) )
            {
                // Cubic I - ibrav=-3
                ibrav = -3;
                // Force to triclinic until special operators coded
                ibrav = 14;
            }
            else
            {
                if ( eqq(abs(a0[2]),abs(a1[2])) && eqq(abs(a1[2]),abs(a2[2])) )
                {
                    // Trigonal 001 axis
                    ibrav =5;
                }
                else
                {
                    // Trigonal, 111 axis
                    ibrav =-5;
                    // Force to triclinic until special operators coded
                    ibrav = 14;
                }
            }
        }
        else if ( eqq(cosab,cosac) && neqq(cosab,cosbc) )
        {
            if ( eqq(abs(a0[0]),abs(a0[1])) && eqq(abs(a1[0]),abs(a1[1])) )
            {
               // Tetragonal I
               ibrav = 7;
               // Force to triclinic until special operators coded
               ibrav = 14;
            }
            else
            {
               // Cubic I - ibrav=3
               ibrav = 3;
            }
        }
        else if ( eqq(cosab,-cosac) && eqq(cosab,cosbc) && eqq(cosab,1.0/3.0) )
        {
            // Cubic I - ibrav=3
            ibrav = 3;
        }
        else if ( eqq(abs(a0[0]),abs(a1[0])) && eqq(abs(a1[1]),abs(a1[1])) )
        {
            // Orthorhombic body-centered
            ibrav = 11;
            // Force to triclinic until special operators coded
            ibrav = 14;
        }
    }
    else if ( eqq(a,b) && neqq(a,c) )
    {
        // Case: a=b/=c
        if ( eqq(cosab,0.0) && eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Case: alpha = beta = gamma = 90
            // Simple tetragonal
            ibrav = 6;
        }
        else if ( eqq(cosab,-0.5) && eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Case: alpha = 120, beta = gamma = 90 => simple hexagonal
            // Simple hexagonal
            ibrav = 4;
        }
        else if ( eqq(cosab,0.5) && eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Case: alpha = 60, beta = gamma = 90 => simple hexagonal
            // Simple hexagonal
            ibrav = 15;
        }
        else if ( eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Orthorhombic bco
            if ( eqq(a0[0],a1[0]) && eqq(a0[1],-a1[1]))
            {
                ibrav = -9;
            }
            else if ( eqq(a0[0],-a1[0]) && eqq(a0[1],a1[1]))
            {
                ibrav = 9;
            }
            // Force to triclinic until special operators coded
            ibrav = 14;
        }
        else if ( eqq(cosac,-cosbc) )
        {
            // bco (unique axis b)
            ibrav =-13;
            // Force to triclinic until special operators coded
            ibrav = 14;
        }
    }
    else if ( eqq(a,c) && neqq(a,b) )
    {
        // Case: a=c/=b
        // Monoclinic bco (unique axis c)
        // QE calls this monoclinic bco but we treat it as Monoclinic
        // ibrav = 13;
        if ( neqq(cosab,0.0) && eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Case: alpha /= 90,  beta = gamma = 90
            // Monoclinic P, unique axis c
            ibrav = 12;
        }
    }
    else if ( eqq(b,c) && neqq(a,b) )
    {
        // Case: a/=b=c
        // Orthorhombic 1-face bco
        if ( eqq(cosab,0.0) && eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Orthorhombic P
            ibrav = 8;
        }
        else
        {
            // Force to triclinic until special operators coded
            ibrav = 14;
        }
    }
    else if ( neqq(a,b) && neqq(a,c) && neqq(b,c) )
    {
        // Case: a/=b/=c
        if ( eqq(cosab,0.0) && eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Case: alpha = beta = gamma = 90
            // Orthorhombic P
            ibrav = 8;
        }
        else if ( neqq(cosab,0.0) && eqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Case: alpha /= 90,  beta = gamma = 90
            // Monoclinic P, unique axis c
            ibrav = 12;
        }
        else if ( eqq(cosab,0.0) && neqq(cosac,0.0) && eqq(cosbc,0.0) )
        {
            // Case: beta /= 90, alpha = gamma = 90
            // Monoclinic P, unique axis b
            ibrav =-12;
            // Force to triclinic until special operators coded
            ibrav = 14;
        }
        else if ( neqq(cosab,0.0) && neqq(cosac,0.0) && neqq(cosbc,0.0) )
        {
            // Case: alpha /= 90, beta /= 90, gamma /= 90
            if ( eqq(abs(a0[0]),abs(a1[0])) && eqq(abs(a0[2]),abs(a2[2])) && eqq(abs(a1[1]),abs(a2[1])) )
            {
                // Orthorhombic F
                ibrav = 10;
                // Force to triclinic until special operators coded
                ibrav = 14;
            }
           else 
           {
              // Triclinic
              ibrav = 14;
           }
        }
    }
    return ibrav;
}

// Saves input vectors in case we rotate them and need the originals
void Lattice::save_vectors(double *a0, double *a1, double *a2)
{
    // Save initial vectors
    for(int i=0;i < 3;i++)a0i[i] = a0[i];
    for(int i=0;i < 3;i++)a1i[i] = a1[i];
    for(int i=0;i < 3;i++)a2i[i] = a2[i];
}

// If the input lattice is one for which RMG has a fast finite diffence representation
// this function will set the lattice vectors so that they coincide with the expected
// RMG orientation for that lattice.
void Lattice::rotate_vectors(double *a0, double *a1, double *a2)
{
    int ibb = lat2ibrav (a0, a1, a2);

    bool dorotate = (ibb == CUBIC_FC) ||
                    (ibb == CUBIC_BC) ||
                    (ibb == ORTHORHOMBIC_PRIMITIVE) ||
                    (ibb == CUBIC_PRIMITIVE) ||
                    (ibb == TETRAGONAL_PRIMITIVE) ||
                    (ibb == HEXAGONAL) ||
                    (ibb == HEXAGONAL2);

    if(!dorotate) return;

    // Get magnitude of vectors
    double m[3] = {0.0, 0.0, 0.0};
    m[0] = sqrt(a0[0]*a0[0] + a0[1]*a0[1] + a0[2]*a0[2]); 
    m[1] = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]); 
    m[2] = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]); 

    for(int i=0;i < 3;i++)a0[i] = 0.0;
    for(int i=0;i < 3;i++)a1[i] = 0.0;
    for(int i=0;i < 3;i++)a2[i] = 0.0;

    double term, cbya;
    switch(ibb)
    {
        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:
            a0[0] = m[0];
            a1[1] = m[1];
            a2[2] = m[2];
            break;
        case CUBIC_FC:
            term = SQRT2 * m[0] / 2.0;
            a0[0] = term;
            a0[1] = term;
            a1[1] = term;
            a1[2] = term;
            a2[0] = term;
            a2[2] = term;
            break;
        case CUBIC_BC:
            term = SQRT3 * m[0] / 2.0;
            for (int ir = 0; ir < 3; ir++)
            {
                a0[ir] = term;
                a1[ir] = term;
                a2[ir] = term;
            }                       /* for ir */
            a0[2] = -term;
            a1[0] = -term;
            a2[1] = -term;
            break;
        case HEXAGONAL:
            cbya = m[2] / m[0];
            a0[0] = m[0];
            a1[0] = -m[0] / 2.0;
            a1[1] =  m[0] * SQRT3 / 2.0;
            a2[2] = m[0] * cbya;
            break;
        case HEXAGONAL2:
            cbya = m[2] / m[0];
            a0[0] = m[0];
            a1[0] = m[0] / 2.0;
            a1[1] =  m[0] * SQRT3 / 2.0;
            a2[2] = m[0] * cbya;
            break;

    }

}


void Lattice::lat2celldm (int ibrav, double alat, double *a1, double *a2, double *a3)
{
/*
       Sets celldm parameters computed from lattice vectors a1,a2,a3 
       a1, a2, a3 are in "alat" units
       If Bravais lattice index ibrav=0, only celldm[0] is set to alat
       See latgen for definition of celldm and lattice vectors.
  
*/
  switch(ibrav)
  {
      case None:
         celldm[0] = 1.0;
         break;
      case CUBIC_PRIMITIVE:
         celldm[0] = sqrt( dot_product (a1,a1) );
         break;
      case CUBIC_FC:
         celldm[0] = sqrt( dot_product (a1,a1) * 2.0 );
         break;
      case -CUBIC_BC:
      case CUBIC_BC:
         celldm[0] = sqrt( dot_product (a1,a1) / 3.0 ) * 2.0;
         break;
      case HEXAGONAL:
         celldm[0] = sqrt( dot_product (a1,a1) );
         celldm[2] = sqrt( dot_product (a3,a3) ) / celldm[0];
         break;
//      case -TRIGONAL_PRIMITIVE:
//      case TRIGONAL_PRIMITIVE:
//         celldm[0] = sqrt( dot_product (a1,a1) );
//         celldm[3] = dot_product(a1,a2) / celldm[0] / sqrt( dot_product( a2,a2) );
//         break;
      case TETRAGONAL_PRIMITIVE:
         celldm[0]= sqrt( dot_product (a1,a1) );
         celldm[2]= sqrt( dot_product (a3,a3) ) / celldm[0];
         break;
//      case TETRAGONAL_BC:
//         celldm[0] = abs(a1(1))*2.0;
//         celldm[2] = abs(a1(3)/a1(1));
//         break;
      case ORTHORHOMBIC_PRIMITIVE:
         celldm[0] = sqrt( dot_product (a1,a1) );
         celldm[1] = sqrt( dot_product (a2,a2) ) / celldm[0];
         celldm[2] = sqrt( dot_product (a3,a3) ) / celldm[0];
         break;
//      case -ORTHORHOMBIC_BASE_CENTRED:
//      case ORTHORHOMBIC_BASE_CENTRED:
//         celldm[0] = abs(a1(1))*2.0;
//         celldm[1] = abs(a2(2))*2.0/celldm[0];
//         celldm[2] = abs(a3(3))/celldm[0];
//         break;
//      case 91:
//         celldm[0] = sqrt( dot_product (a1,a1) );
//         celldm[1] = abs (a2(2))*2.0/celldm[0];
//         celldm[2] = abs (a3(3))*2.0/celldm[0];
//         break;
//      case ORTHORHOMBIC_BC:
//         celldm[0] = abs(a1(1))*2.0;
//         celldm[1] = abs(a2(2))*2.0/celldm[0];
//         celldm[2] = abs(a3(3))*2.0/celldm[0];
//         break;
//      case ORTHORHOMBIC_FC:
//         celldm[0] = abs(a1(1))*2.0;
//         celldm[1] = abs(a1(2))*2.0/celldm[0];
//         celldm[2] = abs(a1(3))*2.0/celldm[0];
//         break;
      case -MONOCLINIC_PRIMITIVE:
      case MONOCLINIC_PRIMITIVE:
         celldm[0] = sqrt( dot_product (a1,a1) );
         celldm[1] = sqrt( dot_product(a2,a2) ) / celldm[0];
         celldm[2] = sqrt( dot_product(a3,a3) ) / celldm[0];
         if ( ibrav == 12 )
         {
            celldm[3] = dot_product(a1,a2) / celldm[0] / sqrt(dot_product(a2,a2));
         }
         else
         {
            celldm[4] = dot_product(a1,a3) / celldm[0] / sqrt(dot_product(a3,a3));
         }
         break;
//      case MONOCLINIC_BASE_CENTRED:
//         celldm[0] = abs(a1(1))*2.0;
//         celldm[1] = sqrt( dot_product(a2,a2)) / celldm[0];
//         celldm[2] = abs (a1(3)/a1(1));
//         celldm[3] = a2(1)/a1(1)/celldm[1]/2.0;
//         break;
//      case -MONOCLINIC_BASE_CENTRED:
//         celldm[0] = abs(a1(1))*2.0;
//         celldm[1] = abs (a2(2)/a2(1));
//         celldm[2] = sqrt( dot_product(a3,a3)) / celldm[0];
//         celldm[4] = a3(1)/a1(1)/celldm[2]/2.0;
//         break;
      case TRICLINIC_PRIMITIVE:
         celldm[0] = sqrt(dot_product(a1,a1));
         celldm[1] = sqrt( dot_product(a2,a2)) / celldm[0];
         celldm[2] = sqrt( dot_product(a3,a3)) / celldm[0];
         celldm[3] = dot_product(a3,a2)/sqrt(dot_product(a2,a2) * dot_product(a3,a3));
         celldm[4] = dot_product(a3,a1) / celldm[0] / sqrt( dot_product(a3,a3));
         celldm[5] = dot_product(a1,a2) / celldm[0] / sqrt(dot_product(a2,a2));
         break;
      default:
         rmg_error_handler (__FILE__, __LINE__, "bravais lattice not programmed.");
      }
      celldm[0] = celldm[0] * alat;
}
