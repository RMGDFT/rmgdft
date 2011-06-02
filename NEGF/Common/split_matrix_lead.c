/*
 **    $Id$    **
 ******************************************************************************/


/*
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>

#include "md.h"
#include "pmo.h"



void split_matrix_lead (int iprobe)
{

    int ntot, i, j, n, ii, jj, li, lj,nstart,idx, index;
    int *desca, ictxt, mb, nb, nprow, npcol, myrow, mycol;
    double distance[9];
    double *H00yz[9], *S00yz[9], *H01yz[9], *S01yz[9], *SCLyz[9], *HCLyz[9];
    double *H00tem, *S00tem, *H01tem, *S01tem, *SCLtem, *HCLtem;
    double blength, clength, yvec, zvec;


    blength = ct.celldm[1] * ct.celldm[0];
    clength = ct.celldm[2] * ct.celldm[0];

    desca = &pmo.desc_lead[(iprobe -1) * DLEN];

    ictxt = desca[1];
    mb = desca[4];
    nb = desca[5];


    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    //  split the matrix for left or right leads, 
    //  left lead is always the same as the first block and the right
    //  lead is always the same as the last block


    if(iprobe == 1) 
    {
        nstart = 0;
        assert(lcr[iprobe].num_states == ct.block_dim[0]);
    }
    else if (iprobe == 2)
    {
        nstart = ct.num_states - lcr[2].num_states;
        assert(lcr[iprobe].num_states == ct.block_dim[ct.num_blocks -1]);
    }
    else
    {
        error_handler(" not programed for other lead in transmission calculation with k points");
        exit(0);
    }    

    H00tem  = lcr[iprobe].H00;
    H01tem  = lcr[iprobe].H01;
    S00tem  = lcr[iprobe].S00;
    S01tem  = lcr[iprobe].S01;
    HCLtem  = lcr[iprobe].HCL;
    SCLtem  = lcr[iprobe].SCL;

    ntot =  pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];

    for(i = 0; i < 9; i++)
    {
        H00yz[i]   = &lcr[iprobe].H00_yz[i * ntot];
        S00yz[i]   = &lcr[iprobe].S00_yz[i * ntot];
        H01yz[i]   = &lcr[iprobe].H01_yz[i * ntot];
        S01yz[i]   = &lcr[iprobe].S01_yz[i * ntot];
        HCLyz[i]   = &lcr[iprobe].HCL_yz[i * ntot];
        SCLyz[i]   = &lcr[iprobe].SCL_yz[i * ntot];
    }        


    for(li = 0; li < pmo.mxllda_lead[iprobe-1]; li++)
    {
        for(lj = 0; lj < pmo.mxlocc_lead[iprobe-1]; lj++)
        {

            /*  li,lj are the  index of distributed matrix */
            /*  i,j are the  index of nondistributed matrix */
            i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
            j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

            /*  ii,jj are the orbital index */
            ii = i+ nstart;
            jj = j+ nstart;

            idx = lj * pmo.mxllda_lead[iprobe-1] + li;

            yvec = states[ii].crds[1] - states[jj].crds[1];
            zvec = states[ii].crds[2] - states[jj].crds[2];

            distance[0] = yvec * yvec + zvec * zvec;
            distance[1] = (yvec + blength) * (yvec + blength) + zvec * zvec;
            distance[2] = (yvec - blength) * (yvec - blength) + zvec * zvec;
            distance[3] = yvec * yvec + (zvec+clength) * (zvec+clength);
            distance[4] = yvec * yvec + (zvec-clength) * (zvec-clength);
            distance[5] = (yvec+blength) * (yvec+blength) + (zvec+clength) * (zvec+clength);
            distance[6] = (yvec+blength) * (yvec+blength) + (zvec-clength) * (zvec-clength);
            distance[7] = (yvec-blength) * (yvec-blength) + (zvec+clength) * (zvec+clength);
            distance[8] = (yvec-blength) * (yvec-blength) + (zvec-clength) * (zvec-clength);

            index =  min_distance_index(distance, 9);

            H00yz[index][idx] = H00tem[idx];
            S00yz[index][idx] = S00tem[idx];
            H01yz[index][idx] = H01tem[idx];
            S01yz[index][idx] = S01tem[idx];
            HCLyz[index][idx] = HCLtem[idx];
            SCLyz[index][idx] = SCLtem[idx];
        }
    }

}
