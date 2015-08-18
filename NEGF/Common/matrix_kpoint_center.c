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


int min_distance_index(double *, int);

void matrix_kpoint_center (complex double *H_tri, double *Stri, double
        *Htri, complex double ene, double kvecy, double kvecz)
{

    int ntot, i, j, n, ii, jj, li, lj,nstart,idx, index;
    int *desca, ictxt, mb, nb, nprow, npcol, myrow, mycol;
    double distance[9], *Htem, *Stem, *Stem1, *Htem1;
    double blength, clength, yvec, zvec;

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

    desca = &pmo.desc_cond[0];

    ictxt = desca[1];
    mb = desca[4];
    nb = desca[5];


    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    //  split the diagonal blocks  
    nstart = 0;
    for(n = 0; n < ct.num_blocks; n++)
    {
        Htem  = &lcr[0].Htri[pmo.diag_begin[n]];
        Stem  = &lcr[0].Stri[pmo.diag_begin[n]];


        for(li = 0; li < pmo.mxllda_cond[n]; li++)
        {
            for(lj = 0; lj < pmo.mxlocc_cond[n]; lj++)
            {

                /*  li,lj are the  index of distributed matrix */
                /*  i,j are the  index of nondistributed matrix */
                i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

                /*  ii,jj are the orbital index */
                ii = i+ nstart;
                jj = j+ nstart;

                idx = lj * pmo.mxllda_cond[n] + li;

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

                H_tri[pmo.diag_begin[n] + idx] = (ene * Stem[idx] - Htem[idx] * Ha_eV) * ctem[index];
            }
        }

        nstart += ct.block_dim[n];

    }



    nstart = 0;
    for(n = 0; n < ct.num_blocks -1; n++)
    {
        Htem  = &lcr[0].Htri[pmo.offdiag_begin[n]];
        Stem  = &lcr[0].Stri[pmo.offdiag_begin[n]];


        for(li = 0; li < pmo.mxllda_cond[n]; li++)
        {
            for(lj = 0; lj < pmo.mxlocc_cond[n+1]; lj++)
            {

                /*  li,lj are the  index of distributed matrix */
                /*  i,j are the  index of nondistributed matrix */
                i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

                /*  ii,jj are the orbital index */
                ii = i+ nstart;
                jj = j+ nstart + ct.block_dim[n];

                idx = lj * pmo.mxllda_cond[n] + li;

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

                H_tri[pmo.offdiag_begin[n] + idx] = (ene * Stem[idx] - Htem[idx] * Ha_eV) * ctem[index];
            }
        }


        nstart += ct.block_dim[n];

    }

//  lower tridigonal block should be ene * conj(S) - conj(H)

    int maxsize, maxrow = 0, maxcol = 0, *descb;
    int ione = 1;
    double one = 1.0, zero = 0.0;

    for(n = 0; n < ct.num_blocks; n++)
    {
        maxrow = rmg_max(maxrow, pmo.mxllda_cond[n]);
        maxcol = rmg_max(maxcol, pmo.mxlocc_cond[n]);
    }

    maxsize = maxrow * maxcol;

    my_malloc_init(Htem1, maxsize, double);
    my_malloc_init(Stem1, maxsize, double);

    nstart = 0;
    for(n = 0; n < ct.num_blocks -1; n++)
    {
        Htem  = &lcr[0].Htri[pmo.offdiag_begin[n]];
        Stem  = &lcr[0].Stri[pmo.offdiag_begin[n]];

        desca = &pmo.desc_cond[( n + (n+1) * ct.num_blocks) * DLEN];
        descb = &pmo.desc_cond[( n+1 + n * ct.num_blocks) * DLEN];
        PDTRAN(&ct.block_dim[n+1], &ct.block_dim[n], &one, Htem, &ione, &ione, desca,
                &zero, Htem1, &ione, &ione, descb);
        PDTRAN(&ct.block_dim[n+1], &ct.block_dim[n], &one, Stem, &ione, &ione, desca,
                &zero, Stem1, &ione, &ione, descb);


        for(li = 0; li < pmo.mxllda_cond[n+1]; li++)
        {
            for(lj = 0; lj < pmo.mxlocc_cond[n]; lj++)
            {

                /*  li,lj are the  index of distributed matrix */
                /*  i,j are the  index of nondistributed matrix */
                i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

                /*  ii,jj are the orbital index */
                ii = i+ nstart + ct.block_dim[n];
                jj = j+ nstart;

                idx = lj * pmo.mxllda_cond[n+1] + li;

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

                H_tri[pmo.lowoffdiag_begin[n] + idx] = (ene * Stem1[idx] - Htem1[idx] * Ha_eV) * ctem[index];
            }
        }


        nstart += ct.block_dim[n];

    }

    my_free(Htem1);
    my_free(Stem1);

}
