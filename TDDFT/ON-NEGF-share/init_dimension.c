/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   main.h for structure defination
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "main.h"

void init_dimension()
{

    int NNBB, NNBBR, NNBBRB;

    NB = 8;
    if (ct.num_states < NB)
        NB = ct.num_states;
    NNBB = (ct.num_states + NB - 1) / NB;
    NNBBR = (NNBB + pct.nprow - 1) / pct.nprow;
    NNBBRB = (NNBBR * NB);
    MXLLDA = NNBBRB ;

    NNBBR = (NNBB + pct.npcol - 1) / pct.npcol;
    MXLCOL = NNBBR * NB;

}


/******/
