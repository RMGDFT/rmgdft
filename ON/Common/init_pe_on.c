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

extern int mpi_nprocs;
extern int mpi_myrank;


void init_pe_on(void)
{

    int npes;
    int ndims, dims[2], periods[2], reorder, kpdelta, remains[2];
    int coords[2];


 /* get total mpi core count  */
    MPI_Comm_size (pct.grid_comm, &npes);

    /* Create a Cartisian topology for parallel in kpoint */
    ndims = 2;
    dims[0] = pct.pe_kpoint;
    dims[1] = npes / dims[0];
    periods[0] = 1;
    periods[1] = 1;
    reorder = 1;
    MPI_Cart_create(pct.img_comm, ndims, &dims[0], &periods[0], reorder, &COMM_KP);

    /* get the coordinate of the processor in COMM_KP */
    MPI_Cart_get(COMM_KP, ndims, &dims[0], &periods[0], &coords[0]);

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
    MPI_Cart_sub(COMM_KP, remains, &COMM_KPSUB1);
    remains[0] = 1;
    remains[1] = 0;
    MPI_Cart_sub(COMM_KP, remains, &COMM_KPSUB2);


    /* XYZ coordinates of this processor */


    my_barrier();
    /* INITIALIZE THE PROCESS GRID */
    int ictxt, info;
    int numst = ct.num_states;
    int rsrc = 0, csrc = 0;
    int mxllda = MXLLDA, nb = ct.scalapack_block_factor;

    sl_init_on(&ictxt, pct.scalapack_nprow, pct.scalapack_npcol);

    Cblacs_gridinfo(ictxt, &pct.scalapack_nprow, &pct.scalapack_npcol, &pct.scalapack_myrow, &pct.scalapack_mycol);



    /* DISTRIBUTE THE MATRIX ON THE PROCESS GRID */
    /* Initialize the array descriptors for the matrices */
    if(pct.scalapack_myrow !=-1)
    {
        DESCINIT(pct.desca, &numst, &numst, &nb, &nb, &rsrc, &csrc, &ictxt, &mxllda, &info);
        if (info != 0)
        {
            printf(" init_pe: DESCINIT, info=%d\n", info);
            fflush(NULL);
            exit(0);

        }
    }

    int numst1;
    numst1 = (ct.num_states + NPES -1)/NPES;
    sl_init_on(&ictxt, 1,NPES);
    int nprow, npcol, myrow,mycol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    DESCINIT(pct.descb, &numst, &numst, &numst, &numst1, &rsrc, &csrc,
            &ictxt, &numst, &info);
    if (info != 0)
    {
        printf(" init_pe for 1xnpes: DESCINIT, info=%d\n", info);
        fflush(NULL);
        exit(0);
    }

    if(myrow ==-1) dprintf("\n WARNNING:  no orbital on processor %d \n", pct.gridpe);
}                               /* end init_pe */


