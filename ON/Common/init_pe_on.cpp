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
#include "transition.h"



extern int mpi_nprocs;
extern int mpi_myrank;


void init_pe_on(void)
{

    MPI_Barrier(pct.img_comm);
    /* INITIALIZE THE PROCESS GRID */
    int ictxt, info;
    int numst = ct.num_states;
    int rsrc = 0, csrc = 0;


    int numst1;
    numst1 = (ct.num_states + pct.grid_npes -1)/pct.grid_npes;
    sl_init_comm(&ictxt, 1, pct.grid_npes, pct.grid_comm);
    int nprow, npcol, myrow,mycol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    descinit(pct.descb, &numst, &numst, &numst, &numst1, &rsrc, &csrc,
            &ictxt, &numst, &info);
    if (info != 0)
    {
        rmg_printf(" init_pe for 1xnpes: descinit, info=%d\n", info);
        fflush(NULL);
        exit(0);
    }

    if(myrow ==-1) dprintf("\n WARNNING:  no orbital on processor %d \n", pct.gridpe);
}                               /* end init_pe */



/*
*  Purpose
*  =======
*
*  SL_INIT initializes an NPROW x NPCOL process grid using a row-major
*  ordering  of  the  processes. This routine retrieves a default system
*  context  which  will  include all available processes. In addition it
*  spawns the processes if needed.
*
*  Arguments
*  =========
*
*  ictxt   (global output) int
*          ictxt specifies the BLACS context handle identifying 
*          created process grid.  The context itself is global.
*
*  NPROW   (global input) int
*          NPROW specifies the number of process rows in the grid
*          to be created.
*
*  NPCOL   (global input) int
*          NPCOL specifies the number of process columns in the grid
*          to be created.
*
*  ============================================================
*  =====================================================================
*/
void sl_init_comm (int *ictxt, int nprow, int npcol, MPI_Comm this_comm)
{
//#if SCALAPACK_LIBS 
    int i, npes;
    int *pmap, *tgmap;

    MPI_Group grp_world, grp_this;

    MPI_Comm_size (this_comm, &npes);
    if (nprow * npcol > npes)
    {
        if(pct.gridpe==0) 
            rmg_printf("Insufficient processes to handle scalapack call, have %d, need %d  * %d", npes, nprow, npcol);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.\n");
    }



    Cblacs_get (0, 0, ictxt);


    /* Allocate space on the assumption that NPES is the same as group size */
    my_malloc (tgmap, npes, int);
    my_malloc (pmap, npes, int);

    /* Set this group rank array maping */
    for (i = 0; i < npes; i++)
        tgmap[i] = i;

    MPI_Comm_group (MPI_COMM_WORLD, &grp_world);
    MPI_Comm_group (this_comm, &grp_this);

    MPI_Group_translate_ranks (grp_this, npes, tgmap, grp_world, pmap);

    Cblacs_gridmap (ictxt, pmap, nprow, nprow, npcol);


    my_free (pmap);
    my_free (tgmap);
//#endif
}

