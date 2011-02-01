/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

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
#include "main.h"
#include "my_scalapack.h"





/*Blocking parameter*/
int nb = 32;


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
void sl_init (int *ictxt, int size)
{
    int i, npes;
    int *pmap, *wgmap, *tgmap;
    int myrow, mycol;
    int nprow, npcol, num_blocks;
    MPI_Group grp_world, grp_this;


    npes = NPES;
    /*First, determine if we want all processors to be involved in scalapack operations
     * We do not want to have too many processors*/
    num_blocks = size / nb;
    if (size % nb)
        num_blocks++;

    /*More processors than 1 per block is a waste */
    if (num_blocks * num_blocks < npes)
        npes = num_blocks * num_blocks;



    /*Distribute processors into 2D grid */
    proc_gridsetup (npes, &nprow, &npcol);

    Cblacs_get (0, 0, ictxt);


    /* calculate MPI world rank range in this group to be mapped to blacs */
    if (nprow * npcol > NPES)
        error_handler
            ("Insufficient MPI group processes to handle scalapack call, have %d, need %d.",
             NPES, nprow * npcol);


    /* Allocate space on the assumption that NPES is the same as group size */
    my_malloc (tgmap, NPES, int);
    my_malloc (pmap, NPES, int);

    /* Set this group rank array maping */
    for (i = 0; (tgmap[i] = i) < NPES; i++);

    /* Get the world rank maping of this groups processes,
       blacs appears to operate on world group ranking */
    MPI_Comm_group (MPI_COMM_WORLD, &grp_world);
    MPI_Comm_group (pct.grid_comm, &grp_this);
    MPI_Group_translate_ranks (grp_this, NPES, tgmap, grp_world, pmap);

    /* Assign nprow*npcol processes to blacs for calculations */
    Cblacs_gridmap (ictxt, pmap, nprow, nprow, npcol);

    /*Figures out blacs information, used so that we can find 
     * which processors are participating in scalapack operations*/
    Cblacs_gridinfo (*ictxt, &nprow, &npcol, &myrow, &mycol);

    /*Variable pct.scalapack_pe will tell use whether the PE participates in scalapack calculations */
    if (myrow >= 0)
        pct.scalapack_pe = 1;
    else
        pct.scalapack_pe = 0;

    /*Store myrow and mycol */
    pct.scalapack_myrow = myrow;
    pct.scalapack_mycol = mycol;

    /*Store number of processor distribution */
    pct.scalapack_nprow = nprow;
    pct.scalapack_npcol = npcol;

    my_free (pmap);

}


void sl_exit (int ictxt)
{
    Cblacs_gridexit (ictxt);
    pct.scalapack_pe = 0;
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


#if !GAMMA_PT
    /*Do not write this info in non-gamma point calculations, the same information is repeated over and over */
    if (!pct.thispe)
        printf ("\n Scalapack processor distribution is nprow:%d and npcol:%d, total processors %d",
                *nprow, *npcol, nproc);
#endif

    if (*npcol * *nprow > nproc)
        error_handler ("\n Wrong processor distribution");

}





/*For given matrix set desca variable, which should be used for matrix operations*/
/*Desca should be declared as      int desca[DLEN]*/
void set_desca (int *desca, int *ictxt, int *size)
{
    int rsrc = 0, csrc = 0, info, izero = 0;
    int mxllda;

    mxllda = NUMROC (size, &nb, &pct.scalapack_myrow, &izero, &pct.scalapack_nprow);
    mxllda = max (1, mxllda);


    if (pct.scalapack_pe)
    {
        /* Initialize the array descriptors for the matrices */
        DESCINIT (desca, size, size, &nb, &nb, &rsrc, &csrc, ictxt, &mxllda, &info);
        if (info)
        {
            printf (" distribute_mat: DESCINIT, info=%d\n", info);
            error_handler ("DESCINIT failed");
        }

    }

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
#if GAMMA_PT
                            dismat[iii + jjj * mxllda] = globmat[ii + jj * size];
#else
                            dismat[2 * (iii + jjj * mxllda)] = globmat[2 * (ii + jj * size)];
                            dismat[2 * (iii + jjj * mxllda) + 1] =
                                globmat[2 * (ii + jj * size) + 1];
#endif
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



#if GAMMA_PT
    for (i = 0; i < size * size; i++)
        globmat[i] = 0.;
#else
    for (i = 0; i < 2 * size * size; i++)
        globmat[i] = 0.;
#endif


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
#if GAMMA_PT
                            globmat[ii + jj * size] = dismat[iii + jjj * mxllda];
#else
                            globmat[2 * (ii + jj * size)] = dismat[2 * (iii + jjj * mxllda)];
                            globmat[2 * (ii + jj * size) + 1] =
                                dismat[2 * (iii + jjj * mxllda) + 1];
#endif
                        }
                    }
                }

            }
        }
    }

}































































#if 0


void diaginit (double *aa, int *desca, double *a, int lda)
{
    int i, ii, iii, li, maxli, iistart;
    int mycol, myrow, nprow, npcol;
    int ictxt = desca[1], mb = desca[4], mxllda = desca[8];



    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    for (i = 0; i < mxllda * mxllda; i++)
        aa[i] = 0.;

    maxli = (lda / (nprow * mb)) + 1;

    if (myrow == mycol)
        for (li = 0; li < maxli; li++)
        {

            iistart = (li * nprow + myrow) * mb;

            for (i = 0; i < mb; i++)
            {

                ii = iistart + i;
                iii = i + li * mb;

                if (iii < mxllda && ii < lda)
                {

                    aa[iii * (mxllda + 1)] = a[ii];
                }

            }
        }

}


/*
*
*     MATGATHER generates and distributes matrice a
*
*/


void distribute_mat (double *bigmat, double *dismat)
{
    int desca[DLEN];
    int ictxt;
    int nb = NB, npcol = NPCOL, nprow = NPROW, numst = ct.num_states;
    int mycol, myrow, mxllda;
    int rsrc = 0, csrc = 0, info;


    mxllda = MXLLDA;
    /* INITIALIZE THE PROCESS GRID */
    sl_init (&ictxt, NPROW, NPCOL);

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);


    /* If I'm in the process grid, execute the program */
    if (myrow != -1)
    {

        /* DISTRIBUTE THE MATRIX ON THE PROCESS GRID */
        /* Initialize the array descriptors for the matrices */
        DESCINIT (desca, &numst, &numst, &nb, &nb, &rsrc, &csrc, &ictxt, &mxllda, &info);
        if (info != 0)
        {
            printf (" distribute_mat: DESCINIT, info=%d\n", info);
            fflush (NULL);
            globalexit (0);
        }


        matinit (dismat, desca, bigmat, ct.num_states);
/*
 *     RELEASE THE PROCESS GRID
 *     Free the BLACS context
 */
        Cblacs_gridexit (ictxt);

    }

}

/********************************************************************/

/********************************************************************/

void get_distributed_mat (double *bigmat, double *dismat)
{
    int desca[DLEN];
    int ictxt;
    int nb = NB, npcol = NPCOL, nprow = NPROW, numst = ct.num_states;
    int mycol, myrow, mxllda;
    int rsrc = 0, csrc = 0, info, idx;
    int n2 = ct.num_states * ct.num_states;

    mxllda = MXLLDA;

    /* INITIALIZE THE PROCESS GRID */
    sl_init (&ictxt, NPROW, NPCOL);

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);


    /* If I'm in the process grid, execute the program */
    if (myrow != -1)
    {

        /* DISTRIBUTE THE MATRIX ON THE PROCESS GRID */
        /* Initialize the array descriptors for the matrices */
        DESCINIT (desca, &numst, &numst, &nb, &nb, &rsrc, &csrc, &ictxt, &mxllda, &info);
        if (info != 0)
        {
            printf (" distribute_mat: DESCINIT, info=%d\n", info);
            fflush (NULL);
            globalexit (0);
        }

        matgather (dismat, desca, bigmat, ct.num_states);
/*
 *     RELEASE THE PROCESS GRID
 *     Free the BLACS context
 */
        Cblacs_gridexit (ictxt);

    }
    else
    {

        for (idx = 0; idx < n2; idx++)
            bigmat[idx] = 0.;

    }


}

/********************************************************************/

void dsymm_dis (char *side, char *uplo, int *nn, double *aa, double *bb, double *cc)
{
    int desca[DLEN];
    int ictxt;
    int nb = NB, npcol = NPCOL, nprow = NPROW;
    int mycol, myrow, mxllda;
    int rsrc = 0, csrc = 0, info;
    _fcd char_fcd1;
    _fcd char_fcd2;
    double zero = 0., one = 1.;
    int ione = 1;

    mxllda = MXLLDA;
    /* INITIALIZE THE PROCESS GRID */
    sl_init (&ictxt, NPROW, NPCOL);

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);


    /* If I'm in the process grid, execute the program */
    if (myrow != -1)
    {

        /* DISTRIBUTE THE MATRIX ON THE PROCESS GRID */
        /* Initialize the array descriptors for the matrices */
        DESCINIT (desca, nn, nn, &nb, &nb, &rsrc, &csrc, &ictxt, &mxllda, &info);
        if (info != 0)
        {
            printf (" distribute_mat: DESCINIT, info=%d\n", info);
            fflush (NULL);
            globalexit (0);
        }
#if CRAY_T3E
        char_fcd1 = _cptofcd (side, 1);
        char_fcd2 = _cptofcd (uplo, 1);
#else
        char_fcd1 = side;
        char_fcd2 = uplo;
#endif

        PSSYMM (char_fcd1, char_fcd2, nn, nn,
                &one, aa, &ione, &ione, desca,
                bb, &ione, &ione, desca, &zero, cc, &ione, &ione, desca);
/*
 *     RELEASE THE PROCESS GRID
 *     Free the BLACS context
 */
        Cblacs_gridexit (ictxt);

    }


}



#if LINUX
void PSSYMM (_fcd side_fcd, _fcd uplo_fcd, int *m, int *n, double *alpha,
             double *a, int *t1, int *t2, int *desca,
             double *b, int *t3, int *t4, int *descb,
             double *beta, double *c, int *t5, int *t6, int *descc)
{
    int lda = ct.num_states, ldb = ct.num_states, ldc = ct.num_states;

    dsymm_ (side_fcd, uplo_fcd, m, n, alpha, a, &lda, b, &ldb, beta, c, &ldc);

}
#endif

#endif
