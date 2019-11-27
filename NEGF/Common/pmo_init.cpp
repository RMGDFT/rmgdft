#include "negf_prototypes.h"
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
#include "init_var.h"
#include "LCR.h"

void pmo_init ()
{

    int ndims, *tgmap, *pmap, dims[2], periods[2], reorder, coords[2], i;
    int myrank, mxllda, *desca, info;
    int iprobe, rsrc =0, csrc =0, numst;
    int idx, nb, mycol, myrow, nprow, npcol;
    int ip, j;
    int izero = 0;
    int remains[2];
    MPI_Group grp_world, grp_this;
   
 
    pmo.npe_energy = pct.grid_npes/pmo.nrow/pmo.ncol; 


    if(pmo.npe_energy * pmo.nrow * pmo.ncol != pct.grid_npes)
    {
        printf("\n parallel matrix grid no good");
        printf("\n pmo.nrow, ncol, pct.grid_npes %d %d %d \n", pmo.nrow, pmo.ncol, pct.grid_npes);
        exit(0);
    }
    pmo.mblock = ct.scalapack_block_factor;


    ndims =2;
    dims[0] = pmo.npe_energy;
    dims[1] = pmo.nrow * pmo.ncol;
    periods[0] = 1;
    periods[1] = 0;
    reorder = 0;

    pmo.orb_index = new int[ct.num_blocks+1];

    my_malloc( pmo.ictxt, pmo.npe_energy, int);
    my_malloc( pmo.mxllda_cond, ct.num_blocks, int);
    my_malloc( pmo.mxlocc_cond, ct.num_blocks, int);
    my_malloc( pmo.diag_begin, ct.num_blocks, int);
    my_malloc( pmo.offdiag_begin, ct.num_blocks, int);
    my_malloc( pmo.lowoffdiag_begin, ct.num_blocks, int);

    my_malloc( pmo.mxllda_lead, cei.num_probe, int);
    my_malloc( pmo.mxlocc_lead, cei.num_probe, int);
    my_malloc( pmo.desc_lead, cei.num_probe * DLEN, int);

    idx = cei.num_probe * ct.num_blocks;
    my_malloc( pmo.desc_cond_lead, idx * DLEN, int);
    my_malloc( pmo.desc_lead_cond, idx * DLEN, int);

    idx = ct.num_blocks * ct.num_blocks;
    my_malloc( pmo.desc_cond, idx * DLEN, int);

    my_malloc(pmap, pmo.nrow * pmo.ncol, int);
    my_malloc(tgmap, pmo.nrow * pmo.ncol, int);

    if( dims[1] * dims[0] != pct.grid_npes) 
    {
        printf("\n processor array not right in pmo_init.c");
        printf("\n npe_energy %d nrow %d, rcol %d\n", pmo.npe_energy, pmo.nrow, pmo.ncol);
        exit(0);
    }


    MPI_Cart_create(pct.grid_comm, ndims, &dims[0], &periods[0], reorder, &COMM_EN);


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

    MPI_Comm_group (MPI_COMM_WORLD, &grp_world);
    MPI_Comm_group (pct.grid_comm, &grp_this);


    for(ip =0; ip < pmo.npe_energy; ip++)
    {
        coords[0] = ip;

        for(i = 0; i < pmo.nrow * pmo.ncol; i++)
        {

            coords[1] = i;
            MPI_Cart_rank( COMM_EN, coords, &myrank);
            tgmap[i] = myrank;
        }


        MPI_Group_translate_ranks (grp_this, dims[1], tgmap, grp_world, pmap);


        Cblacs_get(0, 0, &pmo.ictxt[ip] ); 
        Cblacs_gridmap( &pmo.ictxt[ip], pmap, pmo.nrow, pmo.nrow, pmo.ncol);


    }


    fflush(NULL);
    Cblacs_gridinfo (pmo.ictxt[pmo.myblacs], &nprow, &npcol, &myrow, &mycol);

    fflush(NULL);
    MPI_Barrier(pct.img_comm);
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


        pmo.mxllda_lead[iprobe-1] = numroc(&numst, &pmo.mblock, &myrow, &izero, &nprow);
        pmo.mxlocc_lead[iprobe-1] = numroc(&numst, &pmo.mblock, &mycol, &izero, &npcol);

        mxllda = pmo.mxllda_lead[iprobe -1];

        desca = &pmo.desc_lead[ (iprobe -1) * DLEN];


        /* Initialize the array descriptors for the matrices */
        descinit (desca, &numst, &numst, &nb, &nb, &rsrc, &csrc, &pmo.ictxt[pmo.myblacs], &mxllda, &info);
        if (info != 0)
        {
            printf (" pmo_init: descinit, info=%d\n", info);
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


            mxllda = numroc(&ct.block_dim[i], &pmo.mblock, &myrow, &izero, &nprow);
            desca = &pmo.desc_cond[ idx * DLEN];


            descinit (desca, &ct.block_dim[i], &ct.block_dim[j], &nb, &nb, &rsrc, &csrc, 
                    &pmo.ictxt[pmo.myblacs], &mxllda, &info);
            if (info != 0)
            {
                printf (" distribute_mat: descinit, info=%d\n", info);
                fflush (NULL);
                exit (0);
            } 
        }

    }


    for(i = 0; i < ct.num_blocks; i++)
    {

        pmo.mxllda_cond[i] = numroc(&ct.block_dim[i], &pmo.mblock, &myrow, &izero, &nprow);
        pmo.mxlocc_cond[i] = numroc(&ct.block_dim[i], &pmo.mblock, &mycol, &izero, &npcol);
    }

    pmo.orb_index[0] = 0;

    for (int ib = 0; ib < ct.num_blocks; ib++)
        pmo.orb_index[ib+1] = pmo.orb_index[ib] + ct.block_dim[ib];

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

    pmo.lowoffdiag_begin[0] = pmo.ntot;
    pmo.ntot_low = pmo.ntot + pmo.mxllda_cond[1] * pmo.mxlocc_cond[0];
    for(i = 1; i < ct.num_blocks-1; i++)
    {
        pmo.lowoffdiag_begin[i] = pmo.lowoffdiag_begin[i-1] + pmo.mxllda_cond[i] * pmo.mxlocc_cond[i-1];
        pmo.ntot_low += pmo.mxllda_cond[i+1] * pmo.mxlocc_cond[i];
    }


    for(iprobe = 1; iprobe <=cei.num_probe; iprobe++)
    {
        for(i = 0; i <ct.num_blocks ; i++)
        {
                                                                                                                                            
            idx = i + (iprobe - 1) * ct.num_blocks;
                                                                                                                                            
            mxllda = numroc(&ct.block_dim[i], &pmo.mblock, &myrow, &izero, &nprow);
            numst = lcr[iprobe].num_states;
                                                                                                                                            
            desca = &pmo.desc_cond_lead[ idx * DLEN];
                                                                                                                                            
            descinit (desca, &ct.block_dim[i], &numst, &nb, &nb, &rsrc, &csrc,
                    &pmo.ictxt[pmo.myblacs], &mxllda, &info);
            if (info != 0)
            {
                printf (" distribute_mat: descinit, info=%d\n", info);
                fflush (NULL);
                exit (0);
            }

                                                                                                                                            
            numst = lcr[iprobe].num_states;
            mxllda = numroc(&numst, &pmo.mblock, &myrow, &izero, &nprow);
                                                                                                                                            
            desca = &pmo.desc_lead_cond[ idx * DLEN];
                                                                                                                                            
            descinit (desca, &numst, &ct.block_dim[i], &nb, &nb, &rsrc, &csrc,
                    &pmo.ictxt[pmo.myblacs], &mxllda, &info);
            if (info != 0)
            {
                printf (" distribute_mat: descinit, info=%d\n", info);
                fflush (NULL);
                exit (0);
            }




        }
                                                                                                                                            
    }


/********************************************************/



	my_free(pmap);
} 

