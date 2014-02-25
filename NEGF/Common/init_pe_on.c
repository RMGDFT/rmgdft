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
    ndims = 2;
    dims[0] = pct.pe_kpoint;
    dims[1] = npes / dims[0];
    periods[0] = 1;
    periods[1] = 1;
    reorder = 1;
    MPI_Cart_create (MPI_COMM_WORLD, ndims, &dims[0], &periods[0], reorder, &COMM_KP);

    /* get the coordinate of the processor in COMM_KP */
    MPI_Cart_get (COMM_KP, ndims, &dims[0], &periods[0], &coords[0]);

    pct.coords[0] = coords[0];
    pct.coords[1] = coords[1];
    if (pct.pe_kpoint == 1)
        pct.coords[1] = pct.gridpe;
    /* determine the lower and upper bounds 
     *  of K-point for each group of processors
     */

    kpdelta = (ct.num_kpts + dims[0] - 1) / dims[0];
    pct.kstart = coords[0] * kpdelta;
    pct.kend = pct.kstart + kpdelta;
    if (pct.kend > ct.num_kpts)
        pct.kend = ct.num_kpts;


    /* Partition the communicator into subgroups 
     * COMM_KPSUB1:  dims[0] communicators each with dims[1] processors
     * COMM_KPSUB2:  dims[1] communicators each with dims[0] processors
     */

    remains[0] = 0;
    remains[1] = 1;
    MPI_Cart_sub (COMM_KP, remains, &COMM_KPSUB1);
    remains[0] = 1;
    remains[1] = 0;
    MPI_Cart_sub (COMM_KP, remains, &COMM_KPSUB2);


    /* XYZ coordinates of this processor */

    pe2xyz (pct.coords[1], &ii, &jj, &kk);


    PE_X = pct.pe_x;
    PE_Y = pct.pe_y;
    PE_Z = pct.pe_z;

    /* Now wrap them in case we are running with some processors duplicated */
    /* Two should be enough for any case that we might be doing.            */
    if (ii >= PE_X)
        ii -= PE_X;
    if (ii >= PE_X)
        ii -= PE_X;



    /* Have each processor figure out who it's neighbors are */
    xyz2pe (ii, (jj + 1) % PE_Y, kk, &item);
    coords[1] = item;
    MPI_Cart_rank (COMM_KP, coords, &rank);
    pct.neighbors[NB_N] = rank;

    xyz2pe (ii, (jj - 1 + PE_Y) % PE_Y, kk, &item);
    coords[1] = item;
    MPI_Cart_rank (COMM_KP, coords, &rank);
    pct.neighbors[NB_S] = rank;

    xyz2pe ((ii + 1) % PE_X, jj, kk, &item);
    coords[1] = item;
    MPI_Cart_rank (COMM_KP, coords, &rank);
    pct.neighbors[NB_E] = rank;

    xyz2pe ((ii - 1 + PE_X) % PE_X, jj, kk, &item);
    coords[1] = item;
    MPI_Cart_rank (COMM_KP, coords, &rank);
    pct.neighbors[NB_W] = rank;

    xyz2pe (ii, jj, (kk + 1) % PE_Z, &item);
    coords[1] = item;
    MPI_Cart_rank (COMM_KP, coords, &rank);
    pct.neighbors[NB_U] = rank;

    xyz2pe (ii, jj, (kk - 1 + PE_Z) % PE_Z, &item);
    coords[1] = item;
    MPI_Cart_rank (COMM_KP, coords, &rank);
    pct.neighbors[NB_D] = rank;

    my_barrier ();

    /* Create a Cartisian topology for grid distribution */
    ndims = 3;
    dims[0] = pct.pe_x;
    dims[1] = pct.pe_y;
    dims[2] = pct.pe_z;
    periods[0] = 1;
    periods[1] = 1;
    periods[2] = 1;
    reorder = 0;
    MPI_Cart_create (MPI_COMM_WORLD, ndims, &dims[0], &periods[0], reorder, &COMM_3D);

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
