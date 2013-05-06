/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

/* pmo_init.c initialize the matrix dimensions, BLACK context handle,
 *  and distribution information for every matrix desc...
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "pmo.h"
#include "main.h"

void pmo_init ()
{

    int ndims, *pmap, dims[2], periods[2], reorder, coords[2], *ictxt, i;
    int myrank, mxllda, *desca, info;
    int iprobe, rsrc =0, csrc =0, maxblocks, nblocks, numst;
    int idx, nb, mycol, myrow, nprow, npcol;
    int ip, j;
    int izero = 0;
    int remains[2];
   
 
    pmo.npe_energy = NPES/pmo.nrow/pmo.ncol; 


    if(pmo.npe_energy * pmo.nrow * pmo.ncol != NPES)
    {
        printf("\n parallel matrix grid no good");
        printf("\n pmo.nrow, ncol, NPES %d %d %d \n", pmo.nrow, pmo.ncol, NPES);
        exit(0);
    }
    pmo.mblock = 8;


    ndims =2;
    dims[0] = pmo.npe_energy;
    dims[1] = pmo.nrow * pmo.ncol;
    periods[0] = 1;
    periods[1] = 0;
    reorder = 0;


    my_malloc( pmo.ictxt, pmo.npe_energy, int);
    my_malloc( pmo.mxllda_cond, ct.num_blocks, int);
    my_malloc( pmo.mxlocc_cond, ct.num_blocks, int);
    my_malloc( pmo.diag_begin, ct.num_blocks, int);
    my_malloc( pmo.offdiag_begin, ct.num_blocks, int);

    my_malloc( pmo.mxllda_lead, cei.num_probe, int);
    my_malloc( pmo.mxlocc_lead, cei.num_probe, int);
    my_malloc( pmo.desc_lead, cei.num_probe * DLEN, int);

    idx = cei.num_probe * ct.num_blocks;
    my_malloc( pmo.desc_cond_lead, idx * DLEN, int);
    my_malloc( pmo.desc_lead_cond, idx * DLEN, int);

    idx = ct.num_blocks * ct.num_blocks;
    my_malloc( pmo.desc_cond, idx * DLEN, int);

    my_malloc(pmap, pmo.nrow * pmo.ncol, int);

    if( dims[1] * dims[0] != NPES) 
    {
        printf("\n processor array not right in pmo_init.c");
        printf("\n npe_energy %d nrow %d, rcol %d\n", pmo.npe_energy, pmo.nrow, pmo.ncol);
        exit(0);
    }


    MPI_Cart_create(MPI_COMM_WORLD, ndims, &dims[0], &periods[0], reorder, &COMM_EN);


    MPI_Cart_get(COMM_EN, ndims, dims, periods, coords);


    remains[0] = 1;
    remains[1] = 0;
    MPI_Cart_sub (COMM_EN, remains, &COMM_EN1);
    remains[0] = 0;
    remains[1] = 1;
    MPI_Cart_sub (COMM_EN, remains, &COMM_EN2);

    MPI_Comm_rank(COMM_EN1, &myrank);
    MPI_Comm_rank(COMM_EN2, &myrank);

    pmo.myblacs = coords[0];


    for(ip =0; ip < pmo.npe_energy; ip++)
    {
        coords[0] = ip;

        for(i = 0; i < pmo.nrow * pmo.ncol; i++)
        {

            coords[1] = i;
            MPI_Cart_rank( COMM_EN, coords, &myrank);
            pmap[i] = myrank;
        }

        Cblacs_get(0, 0, &pmo.ictxt[ip] ); 
        Cblacs_gridmap( &pmo.ictxt[ip], pmap, pmo.nrow, pmo.nrow, pmo.ncol);


    }


    fflush(NULL);
    Cblacs_gridinfo (pmo.ictxt[pmo.myblacs], &nprow, &npcol, &myrow, &mycol);

    fflush(NULL);
    my_barrier();
    /* If I'm not in the process grid, return */
    if(myrow == -1) return;



    /* setup the leading dimension for left and right lead */

    for(iprobe = 1; iprobe <=cei.num_probe; iprobe++)
    {

        if( lcr[iprobe].num_states < pmo.mblock )
        {
            pmo.mblock = lcr[iprobe].num_states;
        }
    } 


    for(i = 0; i< ct.num_blocks; i++)
    {
        if(ct.block_dim[i] < pmo.mblock ) 
        {
            pmo.mblock = ct.block_dim[i];
        }
    }


    nb = pmo.mblock;

    /*  initialize for lead matrix */
    for(iprobe = 1; iprobe <=cei.num_probe; iprobe++)
    {

        /* total number of blocks in the row */
        numst = lcr[iprobe].num_states;


        pmo.mxllda_lead[iprobe-1] = NUMROC(&numst, &pmo.mblock, &myrow, &izero, &nprow);
        pmo.mxlocc_lead[iprobe-1] = NUMROC(&numst, &pmo.mblock, &mycol, &izero, &npcol);

        mxllda = pmo.mxllda_lead[iprobe -1];

        desca = &pmo.desc_lead[ (iprobe -1) * DLEN];


        /* Initialize the array descriptors for the matrices */
        DESCINIT (desca, &numst, &numst, &nb, &nb, &rsrc, &csrc, &pmo.ictxt[pmo.myblacs], &mxllda, &info);
        if (info != 0)
        {
            printf (" pmo_init: DESCINIT, info=%d\n", info);
            fflush (NULL);
            exit (0);
        } 

    } 



    /*  initialize for eack blocks in the conductor part */

    /* Initialize the array descriptors for the diagonal block  */
    for(j = 0; j <ct.num_blocks ; j++)
    {
        for(i = 0; i <ct.num_blocks ; i++)
        {

            idx = i + j * ct.num_blocks;


            mxllda = NUMROC(&ct.block_dim[i], &pmo.mblock, &myrow, &izero, &nprow);
            desca = &pmo.desc_cond[ idx * DLEN];

            DESCINIT (desca, &ct.block_dim[i], &ct.block_dim[j], &nb, &nb, &rsrc, &csrc, 
                    &pmo.ictxt[pmo.myblacs], &mxllda, &info);
            if (info != 0)
            {
                printf (" distribute_mat: DESCINIT, info=%d\n", info);
                fflush (NULL);
                exit (0);
            } 
        }

    }


    for(i = 0; i < ct.num_blocks; i++)
    {

        pmo.mxllda_cond[i] = NUMROC(&ct.block_dim[i], &pmo.mblock, &myrow, &izero, &nprow);
        pmo.mxlocc_cond[i] = NUMROC(&ct.block_dim[i], &pmo.mblock, &mycol, &izero, &npcol);
    }


    pmo.diag_begin[0] = 0;

    for(i = 0; i < ct.num_blocks-1; i++)
    {
        pmo.offdiag_begin[i] = pmo.diag_begin[i] + pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
        pmo.diag_begin[i+1] = pmo.offdiag_begin[i] + pmo.mxllda_cond[i] * pmo.mxlocc_cond[i+1];

        if(pct.gridpe == 0) printf (" diag/offdiag begin: %d %d %d \n",
                    i, pmo.offdiag_begin[i], pmo.diag_begin[i+1]); 

    }

    pmo.ntot = pmo.mxllda_cond[0] * pmo.mxlocc_cond[0];
    for(i=1; i< ct.num_blocks; i++)
    {
        pmo.ntot += pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
        pmo.ntot += pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[i];
    }


    /* Initialize the array descriptors for the offdiagonal block  */
	/* which is between the conductor and the lead */
/*
    for(iprobe = 1; iprobe <=cei.num_probe; iprobe++)
    {

        j = cei.probe_in_block[iprobe - 1];

        numst = lcr[iprobe].num_states;
        mxllda = ct.block_dim[j];

        idx = j + (iprobe - 1) * ct.num_blocks;
        desca = &pmo.desc_cond_lead[ idx * DLEN];


    	DESCINIT (desca, &ct.block_dim[j], &numst, &nb, &nb, &rsrc, &csrc, 
                   &pmo.ictxt[pmo.myblacs], &mxllda, &info);
        if (info != 0)
        {
            printf (" distribute_mat: DESCINIT, info=%d\n", info);
            fflush (NULL);
            exit (0);
        } 

    }
*/

    for(iprobe = 1; iprobe <=cei.num_probe; iprobe++)
    {
        for(i = 0; i <ct.num_blocks ; i++)
        {
                                                                                                                                            
            idx = i + (iprobe - 1) * ct.num_blocks;
                                                                                                                                            
            mxllda = NUMROC(&ct.block_dim[i], &pmo.mblock, &myrow, &izero, &nprow);
            numst = lcr[iprobe].num_states;
                                                                                                                                            
            desca = &pmo.desc_cond_lead[ idx * DLEN];
                                                                                                                                            
            DESCINIT (desca, &ct.block_dim[i], &numst, &nb, &nb, &rsrc, &csrc,
                    &pmo.ictxt[pmo.myblacs], &mxllda, &info);
            if (info != 0)
            {
                printf (" distribute_mat: DESCINIT, info=%d\n", info);
                fflush (NULL);
                exit (0);
            }

                                                                                                                                            
            numst = lcr[iprobe].num_states;
            mxllda = NUMROC(&numst, &pmo.mblock, &myrow, &izero, &nprow);
                                                                                                                                            
            desca = &pmo.desc_lead_cond[ idx * DLEN];
                                                                                                                                            
            DESCINIT (desca, &numst, &ct.block_dim[i], &nb, &nb, &rsrc, &csrc,
                    &pmo.ictxt[pmo.myblacs], &mxllda, &info);
            if (info != 0)
            {
                printf (" distribute_mat: DESCINIT, info=%d\n", info);
                fflush (NULL);
                exit (0);
            }




        }
                                                                                                                                            
    }


/********************************************************/



	my_free(pmap);
} 

