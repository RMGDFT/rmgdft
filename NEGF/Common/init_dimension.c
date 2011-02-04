/************************** SVN Revision Information **************************
 **    $Id: init_dimension.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/init_dimension.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 * INPUTS
 * OUTPUT
 *   nothing
 * PARENTS
 *   This is grand-grand-....
 * CHILDREN
 *   run.c
 * SEE ALSO
 *   md.h for structure defination
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "md.h"

void init_dimension ()
{

    int NNBB, NNBBR, NNBBRB;

    NB = 8;
    if (ct.num_states < NB)
        NB = ct.num_states;
    NNBB = (ct.num_states + NB - 1) / NB;
    NNBBR = (NNBB + NPROW - 1) / NPROW;
    NNBBRB = (NNBBR * NB);
    MXLLDA = (NNBBRB + NN % NB);


}


/******/
