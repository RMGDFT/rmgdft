/************************** SVN Revision Information **************************
 **    $Id: distri_fermi.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/distri_fermi.c ************
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 * FUNCTION
 *   void set_energy_weight(real eneR, real eneI, real weight, nenergy)
 *   set up the energies and weights and # of Green functions 
 * INPUTS
 *   
 * OUTPUT
 *   nothing
 * PARENTS
 *   too many
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* This function returns a pointer to a block of memory of size nelem. */
void distri_fermi (REAL eneR, REAL eneI, REAL EF, REAL * distriR, REAL * distriI)
{

    REAL temR, temI, tem1, tem2;

    REAL KT;

    KT = cei.KT;
    temR = (eneR - EF) / KT;
    temI = eneI / KT;
    tem1 = 1.0 + exp (temR) * cos (temI);
    tem2 = exp (temR) * sin (temI);

    *distriR = tem1 / (tem1 * tem1 + tem2 * tem2);
    *distriI = -tem2 / (tem1 * tem1 + tem2 * tem2);
}
