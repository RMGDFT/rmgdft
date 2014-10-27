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
#include "prototypes_on.h"

void init_dimension(int *MXLLDA, int *MXLCOL)
{

    int NNBB, NNBBR, NNBBRB;
    int NB;
    int i, n_factors, *factors, npcol, nprow;

    my_malloc(factors, pct.grid_npes, int);
    n_factors = prime_factors(pct.grid_npes, factors);
    

    nprow = 1;
    npcol = 1;
    for(i = n_factors ; i < 0; i--)
    {
        if(npcol <= nprow) npcol =npcol * factors[i-1];
        if(nprow <  npcol) nprow =nprow * factors[i-1];
    }
    
    if(nprow > npcol)
    {
        pct.scalapack_npcol = nprow;
        pct.scalapack_nprow = npcol;
    }
    else
    {
        pct.scalapack_npcol = npcol;
        pct.scalapack_nprow = nprow;
    }


    my_free(factors);

    
    dprintf("\n scalapcak %d %d\n", pct.scalapack_nprow, pct.scalapack_npcol);
    NB = ct.scalapack_block_factor;
    if (ct.num_states < NB)
        NB = ct.num_states;
    NNBB = (ct.num_states + NB - 1) / NB;
    NNBBR = (NNBB + pct.scalapack_nprow - 1) / pct.scalapack_nprow;
    NNBBRB = (NNBBR * NB);
    *MXLLDA = NNBBRB ;

    NNBBR = (NNBB + pct.scalapack_npcol - 1) / pct.scalapack_npcol;
    *MXLCOL = NNBBR * NB;

}


/******/
