/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
/*
Wenchang Lu, init scalapack ictxt for one communication instead of
usually world 
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"


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
#if SCALAPACK_LIBRARIES
    int i, npes;
    int *pmap, *tgmap;

    MPI_Group grp_world, grp_this;

    MPI_Comm_size (this_comm, &npes);
    if (nprow * npcol > npes)
    {
        error_handler ("Insufficient processes to handle scalapack call, have %d, need %d  * %d", npes, nprow, npcol);
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
#endif
}

