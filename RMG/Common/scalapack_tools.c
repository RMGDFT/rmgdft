/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

/*

	Requires the predefined parameters:
		NN:	size of the global matrix (my_scalapack.h)
		NPES:	number of processors (CFLAGS)

Documentation:

	LAPACK Working Note 94, A User's Guide to the BLACS v1.1

	Manpages: DESCINIT

	WWW: http://www.netlib.org/scalapack/slug/scalapack_slug.html
	(ScaLapack user's guide)



			Fattebert J.-L., November 98 
			fatteber@nemo.physics.ncsu.edu                  

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BaseThread.h"
#include "RmgTimer.h"
#include "grid.h"
#include "main.h"

#if SCALAPACK_LIBS



/*Blocking parameter*/

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
int NB;

static void set_scalapack_comm(int nprow, int npcol, int npes, int images_per_node);
void sl_init (int *ictxt, int size)
{
    int i, npes;
    int *pmap, *tgmap;
    int myrow, mycol;
    int nprow, npcol, num_blocks, izero=0;

    MPI_Group grp_world, grp_this;

    // Allocate the memory for scalapack_desca
    MPI_Comm_size(pct.grid_comm, &npes);
    my_malloc(pct.scalapack_desca, npes*DLEN, int);

    // Reset NB to input value
    NB = ct.scalapack_block_factor;


    /*First, determine if we want all processors to be involved in scalapack operations
     * We do not want to have too many processors*/
    num_blocks = size / NB;
    if (size % NB)
        num_blocks++;

    /*Distribute processors into 2D grid */
    if(npes < ct.images_per_node) 
        error_handler("\n npes %d < images_per_node %d", npes, ct.images_per_node); 
    proc_gridsetup (npes/ct.images_per_node, &nprow, &npcol);

    /*Number of processor in any given direction cannot be more than number of blocks*/
    if(num_blocks < nprow) nprow = num_blocks;
    if(num_blocks < npcol) npcol = num_blocks;


#if GAMMA_PT
    /*Do not write this info in non-gamma point calculations, the same information is repeated over and over */
    printf ("\n Scalapack processor distribution is nprow:%d and npcol:%d, total processors %d",
	    nprow, npcol, nprow*npcol);
#endif
    Cblacs_get (-1, 0, ictxt);

    /* calculate MPI world rank range in this group to be mapped to blacs */
    if (nprow * npcol* ct.images_per_node > npes)
        error_handler
            ("Insufficient MPI group processes to handle scalapack call, have %d, need %d  * %d",
             npes, nprow * npcol, ct.images_per_node);


    /* Allocate space on the assumption that NPES is the same as group size */
    my_malloc (tgmap, npes, int);
    my_malloc (pmap, npes, int);

    /* Set this group rank array maping */
    for (i = 0; i < npes; i++)
	tgmap[i] = i;

    set_scalapack_comm(nprow, npcol, npes, ct.images_per_node);
    /* Get the world rank maping of this groups processes,
       blacs appears to operate on world group ranking */
    int num_pe_scalapack;
    
    MPI_Comm_group (MPI_COMM_WORLD, &grp_world);
    MPI_Comm_group (pct.scalapack_comm, &grp_this);
    MPI_Comm_size (pct.scalapack_comm, &num_pe_scalapack);
    
    MPI_Group_translate_ranks (grp_this, num_pe_scalapack, tgmap, grp_world, pmap);

    /* Assign nprow*npcol processes to blacs for calculations */
    int item; 
    item = pct.thisimg % ct.images_per_node;

    Cblacs_gridmap (ictxt, &pmap[item * nprow * npcol], nprow, nprow, npcol);

    /*Store number of processor distribution */
    pct.scalapack_nprow = nprow;
    pct.scalapack_npcol = npcol;

    /*Figures out blacs information, used so that we can find 
     * which processors are participating in scalapack operations*/
    Cblacs_gridinfo (*ictxt, &nprow, &npcol, &myrow, &mycol);

//dprintf("\n  myrow, mycol nprow npcol %d %d %d %d", myrow, mycol, nprow, npcol);

    /*Variable pct.scalapack_pe will tell use whether the PE participates in scalapack calculations */
    if (myrow >= 0)
        pct.scalapack_pe = 1;
    else
        pct.scalapack_pe = 0;

    /*Store myrow and mycol */
    pct.scalapack_myrow = myrow;
    pct.scalapack_mycol = mycol;



    if(pct.scalapack_pe) {
        pct.scalapack_mpi_rank[myrow*npcol + mycol] = pct.gridpe;
    }
    MPI_Allreduce(MPI_IN_PLACE, pct.scalapack_mpi_rank, npes, MPI_INT, MPI_SUM, pct.grid_comm);

    pct.scalapack_max_dist_size = NUMROC (&ct.num_states, &NB, &myrow, &izero, &nprow) *
                                  NUMROC (&ct.num_states, &NB, &mycol, &izero, &npcol);

    MPI_Allreduce(MPI_IN_PLACE, &pct.scalapack_max_dist_size, 1, MPI_INT, MPI_MAX, pct.grid_comm);

    my_free (pmap);
    my_free (tgmap);

}


void sl_exit (int ictxt, int NotDone)
{
    Cblacs_gridexit (ictxt);
    pct.scalapack_pe = 0;
    Cblacs_exit(NotDone);
    my_free(pct.scalapack_desca);
}




void proc_gridsetup (int nproc, int *nprow, int *npcol)
{

/*This subroutine factorizes the number of processors (nproc)
! into nprow and npcol,  that are the sizes of the 2d processors mesh.
!
! Written by Carlo Cavazzoni
!*/
    int sqrtnp, i;

    sqrtnp = (int) (sqrt (nproc)) + 1;

    for (i = 1; i <= sqrtnp; i++)
        if (nproc % i == 0)
            *nprow = i;

    *npcol = nproc / *nprow;

    if (*npcol * *nprow > nproc)
        error_handler ("\n Wrong processor distribution");

}





/*For given matrix set desca variable, which should be used for matrix operations*/
/*Desca should be declared as      int desca[DLEN]*/
void set_desca (int *desca, int *ictxt, int size)
{
    int rsrc = 0, csrc = 0, info, izero = 0, idx;
    int mxllda, myindex;

    mxllda = NUMROC (&size, &NB, &pct.scalapack_myrow, &izero, &pct.scalapack_nprow);
    mxllda = rmg_max (1, mxllda);


    if (pct.scalapack_pe)
    {
        /* Initialize the array descriptors for the matrices */
        DESCINIT (desca, &size, &size, &NB, &NB, &rsrc, &csrc, ictxt, &mxllda, &info);
        if (info)
        {
            printf (" distribute_mat: DESCINIT, info=%d\n", info);
            error_handler ("DESCINIT failed");
        }

        myindex = pct.scalapack_myrow * pct.scalapack_npcol + pct.scalapack_mycol;
        for(idx = 0;idx < DLEN;idx++) {
            pct.scalapack_desca[myindex*DLEN + idx] = desca[idx];
        }

    }

    int npes;
    MPI_Comm_size(pct.grid_comm, &npes);
    MPI_Allreduce(MPI_IN_PLACE, pct.scalapack_desca, npes*DLEN, MPI_INT, MPI_SUM, pct.grid_comm);
}










/* Sets some stuff and calls function that distributes matrix
 * This requires square matrix*/

void distribute_mat (int *desca, double *bigmat, double *dismat, int *size)
{

    /* If I'm in the process grid, execute the program */
    if (pct.scalapack_pe)
        matinit (desca, dismat, bigmat, *size);
}



/*
*
*     MATINIT generates and distributes matrix
*
*/
void matinit (int *desca, double *dismat, double *globmat, int size)
{
    int i, j, ii, jj, iii, jjj, li, lj, maxrow, maxcol;
    int iistart, jjstart, limb, ljnb, izero = 0;
    int mycol, myrow, nprow, npcol;
    int mb = desca[4], nb = desca[5], mxllda = desca[8];
    int mxlloc = NUMROC (&size, &nb, &pct.scalapack_mycol, &izero, &pct.scalapack_npcol);

    mycol = pct.scalapack_mycol;
    myrow = pct.scalapack_myrow;
    nprow = pct.scalapack_nprow;
    npcol = pct.scalapack_npcol;



    maxrow = (size / (nprow * mb)) + 1;
    maxcol = (size / (npcol * mb)) + 1;


    for (li = 0; li < maxrow; li++)
    {

        iistart = (li * nprow + myrow) * mb;
        limb = li * mb;

        for (lj = 0; lj < maxcol; lj++)
        {

            jjstart = (lj * npcol + mycol) * nb;
            ljnb = lj * nb;

            for (i = 0; i < mb; i++)
            {

                ii = iistart + i;
                iii = i + limb;

                if (iii < mxllda && ii < size)
                {

                    for (j = 0; j < nb; j++)
                    {

                        jj = jjstart + j;
                        jjj = j + ljnb;

                        if (jjj < mxlloc && jj < size)
                        {
                            if(ct.is_gamma) {
                                dismat[iii + jjj * mxllda] = globmat[ii + jj * size];
                            }
                            else {
                                dismat[2 * (iii + jjj * mxllda)] = globmat[2 * (ii + jj * size)];
                                dismat[2 * (iii + jjj * mxllda) + 1] =
                                    globmat[2 * (ii + jj * size) + 1];
                            }
                        }
                    }
                }

            }
        }
    }

}






void matgather (double *dismat, int *desca, double *globmat, int size)
{
    int i, j, ii, jj, iii, jjj, li, lj, maxcol, maxrow, jjstart, iistart;
    int limb, ljnb, izero = 0;
    int mycol, myrow, nprow, npcol;
    int mb = desca[4], nb = desca[5], mxllda = desca[8];
    int mxlloc = NUMROC (&size, &nb, &pct.scalapack_mycol, &izero, &pct.scalapack_npcol);


    mycol = pct.scalapack_mycol;
    myrow = pct.scalapack_myrow;
    nprow = pct.scalapack_nprow;
    npcol = pct.scalapack_npcol;



    if(ct.is_gamma) {
        for (i = 0; i < size * size; i++)
            globmat[i] = 0.;
    }
    else {
        for (i = 0; i < 2 * size * size; i++)
            globmat[i] = 0.;
    }


    maxrow = (size / (nprow * mb)) + 1;
    maxcol = (size / (npcol * mb)) + 1;

    /* loop on the blocks of size mb*nb */
    for (li = 0; li < maxrow; li++)
    {

        iistart = (li * nprow + myrow) * mb;
        limb = li * mb;

        for (lj = 0; lj < maxcol; lj++)
        {

            jjstart = (lj * npcol + mycol) * nb;
            ljnb = lj * nb;

            /* loop in the block */
            for (i = 0; i < mb; i++)
            {

                ii = iistart + i;
                iii = limb + i;

                if (iii < mxllda && ii < size)
                {

                    for (j = 0; j < nb; j++)
                    {

                        jj = jjstart + j;
                        jjj = j + ljnb;

                        //if (jjj < mxllda && jj < size)
                        if (jjj < mxlloc && jj < size)
                        {
                            if(ct.is_gamma) {
                                globmat[ii + jj * size] = dismat[iii + jjj * mxllda];
                            }
                            else {
                                globmat[2 * (ii + jj * size)] = dismat[2 * (iii + jjj * mxllda)];
                                globmat[2 * (ii + jj * size) + 1] =
                                    dismat[2 * (iii + jjj * mxllda) + 1];
                            }
                        }
                    }
                }

            }
        }
    }

}


// Used to sum a distributed matrix over all nodes and store the results
// to the scalapack PE's.
//
void matsum (double *tmat, double *dismat, double *globmat, int size)
{
    int idx, jdx;
    int istop, len;
    int col, row, nprow, npcol, dist_length[MAX_RMG_THREADS];

    nprow = pct.scalapack_nprow;
    npcol = pct.scalapack_npcol;

    istop = (nprow*npcol) / ct.THREADS_PER_NODE;
    istop = istop * ct.THREADS_PER_NODE;
    for(idx = 0;idx < istop;idx+=ct.THREADS_PER_NODE) {

// Use OpenMP for the packing within a node
#pragma omp parallel for private(row,col)
        for(jdx = 0;jdx < ct.THREADS_PER_NODE;jdx++) {
            row = (idx + jdx) / npcol;
            col = (idx + jdx) % npcol;
            dist_length[jdx] = matsum_packbuffer(row, col, &tmat[jdx * pct.scalapack_max_dist_size], globmat, size);
        }

        for(jdx = 0;jdx < ct.THREADS_PER_NODE;jdx++) {
            row = (idx + jdx) / npcol;
            col = (idx + jdx) % npcol;
            // Reduce to the target
            MPI_Reduce(&tmat[jdx * pct.scalapack_max_dist_size], dismat, dist_length[jdx], MPI_DOUBLE, 
                       MPI_SUM, pct.scalapack_mpi_rank[row*npcol + col], pct.grid_comm);
        }

    }

    // Finish up the remainder
    for(idx = istop;idx < nprow*npcol;idx++) {
        row = idx / npcol;
        col = idx % npcol;
        len = matsum_packbuffer(row, col, tmat, globmat, size);
        MPI_Reduce(tmat, dismat, len, MPI_DOUBLE, MPI_SUM, pct.scalapack_mpi_rank[row*npcol + col], pct.grid_comm);
    }

#if 0
    for(row =0;row < nprow;row++) {

        for(col =0;col < npcol;col++) {

            len = matsum_packbuffer(row, col, tmat, globmat, size);

            // Reduce to the target node
            MPI_Reduce(tmat, dismat, len, MPI_DOUBLE, MPI_SUM, pct.scalapack_mpi_rank[row*npcol + col], pct.grid_comm);

        }

    }
#endif

}

int matsum_packbuffer(int row, int col, double *buffer, double *globmat, int size)
{

    int i, j, ii, jj, iii, jjj, li, lj, maxrow, maxcol;
    int iistart, jjstart, limb, ljnb, izero = 0;
    int nprow, npcol, dist_length;
    int mb, nb, mxllda, mxlloc;

    nprow = pct.scalapack_nprow;
    npcol = pct.scalapack_npcol;
    mb = pct.scalapack_desca[(row*npcol + col)*DLEN + 4];
    nb = pct.scalapack_desca[(row*npcol + col)*DLEN + 5];
    maxrow = (size / (nprow * mb)) + 1;
    maxcol = (size / (npcol * mb)) + 1;
    mxllda = pct.scalapack_desca[(row*npcol + col)*DLEN + 8];
    mxlloc = NUMROC (&size, &nb, &col, &izero, &pct.scalapack_npcol);
    dist_length = NUMROC (&ct.num_states, &mb, &row, &izero, &pct.scalapack_nprow) * 
                  NUMROC (&ct.num_states, &mb, &col, &izero, &pct.scalapack_npcol);
#if !GAMMA_PT
    dist_length *= 2;
#endif

    for (li = 0; li < maxrow; li++)
    {

        iistart = (li * nprow + row) * mb;
        limb = li * mb;

        for (lj = 0; lj < maxcol; lj++)
        {

            jjstart = (lj * npcol + col) * nb;
            ljnb = lj * nb;

            for (i = 0; i < mb; i++)
            {

                ii = iistart + i;
                iii = i + limb;

                if (iii < mxllda && ii < size)
                {

                    for (j = 0; j < nb; j++)
                    {

                        jj = jjstart + j;
                        jjj = j + ljnb;

                        if (jjj < mxlloc && jj < size)
                        {


#if GAMMA_PT
                            buffer[iii + jjj * mxllda] = globmat[ii + jj * size];
#else
                            buffer[2 * (iii + jjj * mxllda)] = globmat[2 * (ii + jj * size)];
                            buffer[2 * (iii + jjj * mxllda) + 1] =
                                globmat[2 * (ii + jj * size) + 1];
#endif

                        }
                    }
                }

            }
        }
    }

    return dist_length;

}


// Performs a reduce and distribute operation to the scalapack PE's
// global_matrix  contains the data to be reduced while dist_matrix
// is the target distributed matrix and work is temporary space that
// should be sized at least as large as dist_matrix
void reduce_and_dist_matrix(int n, double *global_matrix, double *dist_matrix, double *work)
{

    int stop;
    void *RT = BeginRmgTimer("Diagonalization: reduce and distribute");
    stop = n * n;

    if(ct.scalapack_global_sums) {

        /*Sum matrix over all processors */
        void *RT1 = BeginRmgTimer("Diagonalization: MPI_Allreduce");
        MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        EndRmgTimer(RT1);

        /*Distribute global matrix A */
        if (pct.scalapack_pe)
            distribute_mat (pct.desca, global_matrix, dist_matrix, &n);

    }
    else {

        matsum (work, dist_matrix, global_matrix, n);

    }

    EndRmgTimer(RT);
}

static void set_scalapack_comm(int nprow, int npcol, int npes, int images_per_node)
{

    //devide NPES to groups, each group has same number of MPI processes which is >= nprw * npcol * imagpes_per_nodes.
    int num_group, i, num_pe;
    MPI_Comm tem_comm;

    num_group = npes /(nprow * npcol * images_per_node);

    for ( i = num_group; i >0; i--)
        if(npes % i == 0) break;

    num_group = i;
    num_pe = npes/i;
    
    pct.scalapack_npes = num_pe;

    int ndims = 2;
    int dims[] = { num_pe, num_group};
    int periods[] = { 0, 0 };
    int reorder = 1;
    int remains[] = { 1, 0 };

    MPI_Cart_create (pct.grid_comm, ndims, dims, periods, reorder, &tem_comm);
    MPI_Cart_sub (tem_comm, remains, &pct.scalapack_comm);


}
#endif
