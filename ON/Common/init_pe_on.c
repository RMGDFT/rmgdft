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



#include "md.h"
#include <stdlib.h>
#include <stdio.h>


extern int mpi_nprocs;
extern int mpi_myrank;


void init_pe_on(void)
{

    int npes, ii, jj, kk;
    int ndims, dims[2], periods[2], reorder, kpdelta, remains[2];
    int coords[2], item, rank, PE_X, PE_Y, PE_Z;


 /* get total mpi core count  */
    MPI_Comm_size (pct.grid_comm, &npes);

    /* Create a Cartisian topology for parallel in kpoint */
    ndims = 2;
    dims[0] = pct.pe_kpoint;
    dims[1] = npes / dims[0];
    periods[0] = 1;
    periods[1] = 1;
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, &dims[0], &periods[0], reorder, &COMM_KP);

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

    pe2xyz(pct.coords[1], &ii, &jj, &kk);


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
    xyz2pe(ii, (jj + 1) % PE_Y, kk, &item);
    coords[1] = item;
    MPI_Cart_rank(COMM_KP, coords, &rank);
    pct.neighbors[NB_N] = rank;

    xyz2pe(ii, (jj - 1 + PE_Y) % PE_Y, kk, &item);
    coords[1] = item;
    MPI_Cart_rank(COMM_KP, coords, &rank);
    pct.neighbors[NB_S] = rank;

    xyz2pe((ii + 1) % PE_X, jj, kk, &item);
    coords[1] = item;
    MPI_Cart_rank(COMM_KP, coords, &rank);
    pct.neighbors[NB_E] = rank;

    xyz2pe((ii - 1 + PE_X) % PE_X, jj, kk, &item);
    coords[1] = item;
    MPI_Cart_rank(COMM_KP, coords, &rank);
    pct.neighbors[NB_W] = rank;

    xyz2pe(ii, jj, (kk + 1) % PE_Z, &item);
    coords[1] = item;
    MPI_Cart_rank(COMM_KP, coords, &rank);
    pct.neighbors[NB_U] = rank;

    xyz2pe(ii, jj, (kk - 1 + PE_Z) % PE_Z, &item);
    coords[1] = item;
    MPI_Cart_rank(COMM_KP, coords, &rank);
    pct.neighbors[NB_D] = rank;

    my_barrier();
    /* INITIALIZE THE PROCESS GRID */
    int ictxt, info;
    int numst = ct.num_states;
    int rsrc = 0, csrc = 0;
    int mxllda = MXLLDA, nb = NB;

    sl_init(&ictxt, pct.nprow, pct.npcol);

    Cblacs_gridinfo(ictxt, &pct.nprow, &pct.npcol, &pct.myrow, &pct.mycol);



    /* DISTRIBUTE THE MATRIX ON THE PROCESS GRID */
    /* Initialize the array descriptors for the matrices */
    if(pct.myrow !=-1)
    {
        DESCINIT(pct.desca, &numst, &numst, &nb, &nb, &rsrc, &csrc, &ictxt, &mxllda, &info);
        if (info != 0)
        {
            printf(" init_pe: DESCINIT, info=%d\n", info);
            fflush(NULL);
            exit(0);

        }
    }

    int numst1, npes_tem;
    numst1 = (ct.num_states + NPES -1)/NPES;
    npes_tem = (ct.num_states + numst1 -1) /numst1;
    sl_init(&ictxt, 1,npes_tem);
    int nprow, npcol, myrow,mycol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    DESCINIT(pct.descb, &numst, &numst, &numst, &numst1, &rsrc, &csrc,
            &ictxt, &numst, &info);

}                               /* end init_pe */


