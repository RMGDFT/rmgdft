/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/init_pe.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void init_pe(void)
 *   Initializes processor control structures.
 *   Make each processor knows who are neighbors
 * INPUTS
 *   nothing
 * OUTPUT
 *   neighbor informations are stored in pct.neighbors
 * PARENTS
 *   md.c
 * CHILDREN
 *   pe2xyz.c  
 * SOURCE
 */



#include "main.h"
#include "prototypes_on.h"
#include <stdlib.h>
#include <stdio.h>
#include "init_var.h"
#include "blacs.h"


extern int mpi_nprocs;
extern int mpi_myrank;


void init_pe_on(void)
{

    MPI_Barrier(pct.img_comm);
    /* INITIALIZE THE PROCESS GRID */
    int ictxt, info;
    int numst = ct.num_states;
    int rsrc = 0, csrc = 0;
    int mxllda = MXLLDA, nb = ct.scalapack_block_factor;

    sl_init_comm(&ictxt, pct.scalapack_nprow, pct.scalapack_npcol, pct.grid_comm);

    Cblacs_gridinfo(ictxt, &pct.scalapack_nprow, &pct.scalapack_npcol, &pct.scalapack_myrow, &pct.scalapack_mycol);



    /* DISTRIBUTE THE MATRIX ON THE PROCESS GRID */
    /* Initialize the array descriptors for the matrices */
    if(pct.scalapack_myrow !=-1)
    {
        descinit(pct.desca, &numst, &numst, &nb, &nb, &rsrc, &csrc, &ictxt, &mxllda, &info);
        if (info != 0)
        {
            printf(" init_pe: descinit, info=%d\n", info);
            fflush(NULL);
            exit(0);

        }
    }

    int numst1;
    numst1 = (ct.num_states + pct.grid_npes -1)/pct.grid_npes;
    sl_init_comm(&ictxt, 1, pct.grid_npes, pct.grid_comm);
    int nprow, npcol, myrow,mycol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    descinit(pct.descb, &numst, &numst, &numst, &numst1, &rsrc, &csrc,
            &ictxt, &numst, &info);
    if (info != 0)
    {
        printf(" init_pe for 1xnpes: descinit, info=%d\n", info);
        fflush(NULL);
        exit(0);
    }

    if(myrow ==-1) dprintf("\n WARNNING:  no orbital on processor %d \n", pct.gridpe);
}                               /* end init_pe */


