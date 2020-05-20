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
#include "Scalapack.h"

void init_dimension(int *MXLLDA, int *MXLCOL)
{

    Scalapack *MainSp;
    int scalapack_groups = 1;
    int numst = ct.num_states;
    int last = 1;
    MainSp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, numst,
            ct.scalapack_block_factor, last, pct.grid_comm);

    
    pct.scalapack_nprow = MainSp->GetRows();
    pct.scalapack_npcol = MainSp->GetCols();
    pct.scalapack_myrow = MainSp->GetRow();
    pct.scalapack_mycol = MainSp->GetCol();
    
    int *desca = MainSp->GetDistDesca();
    
    for(int i = 0; i < DLEN; i++) pct.desca[i] = desca[i];

    *MXLLDA = MainSp->GetDistMdim();

    *MXLCOL = MainSp->GetDistNdim();

}


/******/
