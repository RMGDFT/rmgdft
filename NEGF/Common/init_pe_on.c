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
#include "init_var.h"
#include "LCR.h"
#include <stdlib.h>
#include <stdio.h>


extern int mpi_nprocs;
extern int mpi_myrank;


void init_pe_on (void)
{

    int npes, ii, jj, kk;
    int ndims, dims[3], periods[3], reorder, kpdelta, remains[3];
    int coords[3], item, rank, PE_X, PE_Y, PE_Z;


 /* get total mpi core count  */
    MPI_Comm_size (pct.grid_comm, &npes);

    /* Create a Cartisian topology for parallel in kpoint */



    /* Create a Cartisian topology for grid distribution */
    ndims = 3;
    dims[0] = pct.pe_x;
    dims[1] = pct.pe_y;
    dims[2] = pct.pe_z;
    periods[0] = 1;
    periods[1] = 1;
    periods[2] = 1;
    reorder = 0;
    MPI_Cart_create (pct.grid_comm, ndims, &dims[0], &periods[0], reorder, &COMM_3D);

    /* Partition the communicator COMM_3D  into subgroups 
     * COMM_PEX:  pey * pez communicators each with pex processors
     * COMM_PEY:  pex * pez communicators each with pey processors
     * COMM_PEZ:  pey * pex communicators each with pez processors
     */

    remains[0] = 1;
    remains[1] = 0;
    remains[2] = 0;
    MPI_Cart_sub (COMM_3D, remains, &COMM_PEX);

    remains[0] = 0;
    remains[1] = 1;
    remains[2] = 0;
    MPI_Cart_sub (COMM_3D, remains, &COMM_PEY);

    remains[0] = 0;
    remains[1] = 0;
    remains[2] = 1;
    MPI_Cart_sub (COMM_3D, remains, &COMM_PEZ);

    pe2xyz (pct.gridpe, &ii, &jj, &kk);
    MPI_Comm_size (COMM_PEX, &PE_Y);
    MPI_Comm_rank (COMM_PEX, &PE_X);

    MPI_Cart_get (COMM_3D, ndims, &dims[0], &periods[0], &coords[0]);

}                               /* end init_pe */


/******/
