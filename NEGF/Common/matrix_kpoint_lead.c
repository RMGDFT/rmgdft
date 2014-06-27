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

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"



void matrix_kpoint_lead (complex double *S00, complex double *H00,
        complex double *S01, complex double *H01, complex double *SCL, complex
        double *HCL, double kvecy, double kvecz, int iprobe)
{

    int ntot, i, j, n, ii, jj, li, lj,nstart,idx, index;
    int *desca, ictxt, mb, nb, nprow, npcol, myrow, mycol;
    double distance[9];
    double blength, clength, yvec, zvec;

    int nstart_block, subsystem, idx0, llda, locc;

    double complex  ctem[9];

    ctem[0] = 1.0;
    ctem[1] = cexp(+I * kvecy);
    ctem[2] = cexp(-I * kvecy);
    ctem[3] = cexp(+I * kvecz);
    ctem[4] = cexp(-I * kvecz);
    ctem[5] = cexp(+I * kvecy + I * kvecz);
    ctem[6] = cexp(+I * kvecy - I * kvecz);
    ctem[7] = cexp(-I * kvecy + I * kvecz);
    ctem[8] = cexp(-I * kvecy - I * kvecz);


    blength = get_celldm(1) * get_celldm(0);
    clength = get_celldm(2) * get_celldm(0);

    desca = &pmo.desc_lead[(iprobe -1) * DLEN];

    ictxt = desca[1];
    mb = desca[4];
    nb = desca[5];


    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    //  split the matrix for left or right leads, 
    //  left lead is always the same as the first block and the right
    //  lead is always the same as the last block

    nstart = 0;
    for (subsystem = 0; subsystem < cei.num_subsystem; subsystem++)
    {
        idx0 = cei.subsystem_idx[subsystem];
        if(idx0 == iprobe) break;
        nstart += lcr[idx0].state_end - lcr[idx0].state_begin;
    }



    llda =  pmo.mxllda_lead[iprobe-1] ;
    locc =  pmo.mxlocc_lead[iprobe-1];
    ntot =  llda * locc;

    for(i = 0; i < ntot; i++) 
    {
        S01[i] = 0.0;
        H01[i] = 0.0;
        S00[i] = 0.0;
        H00[i] = 0.0;
        SCL[i] = 0.0;
        HCL[i] = 0.0;
    }

    for(li = 0; li < llda; li++)
    {
        for(lj = 0; lj < locc; lj++)
        {

            /*  li,lj are the  index of distributed matrix */
            /*  i,j are the  index of nondistributed matrix */
            i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
            j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

            /*  ii,jj are the orbital index */
            ii = i+ nstart;
            jj = j+ nstart;

            idx = lj * llda + li;

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

            S01[idx] += lcr[iprobe].S01[idx] * ctem[index];
            H01[idx] += lcr[iprobe].H01[idx] * ctem[index];
            S00[idx] += lcr[iprobe].S00[idx] * ctem[index];
            H00[idx] += lcr[iprobe].H00[idx] * ctem[index];

        }
    }


    i = cei.probe_in_block[iprobe - 1];

    llda  = pmo.mxllda_cond[i];
    locc  = pmo.mxlocc_lead[iprobe-1];
    ntot  = llda * locc;

    nstart_block = 0;
    for (j =0; j < cei.probe_in_block[iprobe-1]; j++)
        nstart_block += ct.block_dim[j];


    //    printf("\n iprobe nstart %d %d %d ", iprobe, nstart, nstart_block);


    for(li = 0; li < llda; li++)
    {
        for(lj = 0; lj < locc; lj++)
        {

            /*  li,lj are the  index of distributed matrix */
            /*  i,j are the  index of nondistributed matrix */
            i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
            j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

            /*  ii,jj are the orbital index */
            ii = i+ nstart_block;
            jj = j+ nstart;

            idx = lj * llda + li;

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

            SCL[idx] += lcr[iprobe].SCL[idx] * ctem[index];
            HCL[idx] += lcr[iprobe].HCL[idx] * ctem[index];
        }
    }


}
