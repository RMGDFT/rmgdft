/************************** SVN Revision Information **************************
 **    $Id: tools_submat.c 1242 2011-02-02 18:55:23Z luw $    **
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
#include <float.h>
#include <assert.h>
#include "md.h"
#include "my_scalapack.h"

#if USE_DIS_MAT


#define globalexit  exit


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
void sl_init (int *ictxt, int nprow, int npcol)
{
    int i, iam, nprocs;
    /*char    order='R'; */
    int *pmap;

    assert (nprow > 0);
    assert (npcol > 0);

    Cblacs_pinfo (&iam, &nprocs);

/*
 *  If machine needs additional set up, do it now
 */

/*  Wenchang LU
    if( nprocs<=1 ){
        if( iam==0)nprocs=nprow*npcol;
        Cblacs_setup(&iam, &nprocs);
    }
*/
/*
 *  Define process grid
 */
/*   Wenchang Lu
    Cblacs_get(-1, 0, ictxt);
*/
    Cblacs_get (0, 0, ictxt);

    /*Cblacs_gridinit(ictxt, &order, nprow, npcol); */
    /*
       iam=0;
       for(i=0;i<nprow;i++)
       for(j=0;j<npcol;j++){
       pmap[i][j]=iam;
       iam++;
       }
     */


    my_malloc( pmap, NPROW * NPCOL, int );
    for (i = 0; i < nprow * npcol; i++)
        pmap[i] = i;

    Cblacs_gridmap (ictxt, pmap, nprow, nprow, npcol);
    my_free(pmap);

}


/*
*
*     MATINIT generates and distributes matrice a
*
*/
void matinit (double *aa, int *desca, double *a, int lda)
{
    int i, j, ii, jj, iii, jjj, li, lj, maxli;
    int iistart, jjstart, limb, ljnb;
    int mycol, myrow, nprow, npcol;
    int ictxt = desca[1], mb = desca[4], nb = desca[5], mxllda = desca[8];



    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);


    maxli = (lda / (nprow * mb)) + 1;


    for (li = 0; li < maxli; li++)
    {

        iistart = (li * nprow + myrow) * mb;
        limb = li * mb;

        for (lj = 0; lj < maxli; lj++)
        {

            jjstart = (lj * npcol + mycol) * nb;
            ljnb = lj * nb;

            for (i = 0; i < mb; i++)
            {

                ii = iistart + i;
                iii = i + limb;

                if (iii < mxllda && ii < lda)
                {

                    for (j = 0; j < nb; j++)
                    {

                        jj = jjstart + j;
                        jjj = j + ljnb;

                        if (jjj < mxllda && jj < lda)
                            aa[iii + jjj * mxllda] = a[ii + jj * lda];
                    }
                }

            }
        }
    }

}


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
void matgather (double *aa, int *desca, double *a, int lda)
{
    int i, j, ii, jj, iii, jjj, li, lj, maxli, jjstart, iistart;
    int limb, ljnb;
    int mycol, myrow, nprow, npcol;
    int ictxt = desca[1], mb = desca[4], nb = desca[5], mxllda = desca[8];



    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    for (i = 0; i < lda * lda; i++)
        a[i] = 0.;


    maxli = (lda / (nprow * mb)) + 1;

    /* loop on the blocks of size mb*nb */
    for (li = 0; li < maxli; li++)
    {

        iistart = (li * nprow + myrow) * mb;
        limb = li * mb;

        for (lj = 0; lj < maxli; lj++)
        {

            jjstart = (lj * npcol + mycol) * nb;
            ljnb = lj * nb;

            /* loop in the block */
            for (i = 0; i < mb; i++)
            {

                ii = iistart + i;
                iii = limb + i;

                if (iii < mxllda && ii < lda)
                {

                    for (j = 0; j < nb; j++)
                    {

                        jj = jjstart + j;
                        jjj = j + ljnb;

                        if (jjj < mxllda && jj < lda)
                            a[ii + jj * lda] = aa[iii + jjj * mxllda];
                    }
                }

            }
        }
    }

}



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
void print_distribute_mat (double *dismat)
{
    int desca[DLEN];
    int ictxt;
    int nb = NB, npcol = NPCOL, nprow = NPROW, numst = ct.num_states;
    int mycol, myrow, mxllda;
    int rsrc = 0, csrc = 0, info;
    int n2 = ct.num_states * ct.num_states, idx;


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

    }
    for (idx = 0; idx < n2; idx++)
        work_matrix[idx] = 0.;
    if (myrow != -1)
        matgather (dismat, desca, work_matrix, numst);
    global_sums (work_matrix, &n2);
    if (pct.thispe == 0)
    {
        printf (" Distributed matrix\n");
        print_matrix (work_matrix, 5, ct.num_states);
    }
    if (myrow != -1)
    {
/*
 *     RELEASE THE PROCESS GRID
 *     Free the BLACS context
 */
        Cblacs_gridexit (ictxt);

    }

    fflush (NULL);
    my_barrier ();


}

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

#endif /* USE_DIS_MAT */


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
