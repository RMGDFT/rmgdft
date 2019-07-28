#include "negf_prototypes.h"
/*
 **    $Id: matrix_kpoint_center.c 2381 2014-06-30 19:35:11Z luw $    **
******************************************************************************/
 

/*
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"
#include "blacs.h"


int min_distance_index(double *, int);

void green_kpoint_phase (std::complex<double> *green, double kvecy, double kvecz, int up_and_low)
{

    int ntot, i, j, n, ii, jj, li, lj,nstart,idx, index;
    int *desca, ictxt, mb, nb, nprow, npcol, myrow, mycol;
    double distance[9], *Htem, *Stem, *Stem1, *Htem1;
    double blength, clength, yvec, zvec;

   std::complex<double>  ctem[9];
   std::complex<double> I(0.0, 1.0);

    ctem[0] = 1.0;
    ctem[1] = std::exp(+I * kvecy);
    ctem[2] = std::exp(-I * kvecy);
    ctem[3] = std::exp(+I * kvecz);
    ctem[4] = std::exp(-I * kvecz);
    ctem[5] = std::exp(+I * kvecy + I * kvecz);
    ctem[6] = std::exp(+I * kvecy - I * kvecz);
    ctem[7] = std::exp(-I * kvecy + I * kvecz);
    ctem[8] = std::exp(-I * kvecy - I * kvecz);

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

                green[pmo.diag_begin[n] + idx] *= conj(ctem[index]);
            }
        }

        nstart += ct.block_dim[n];

    }



    nstart = 0;
    for(n = 0; n < ct.num_blocks -1; n++)
    {


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

                green[pmo.offdiag_begin[n] + idx] *= conj(ctem[index]);
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

    if(up_and_low == 1)
    {
        nstart = 0;
        for(n = 0; n < ct.num_blocks -1; n++)
        {

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

                    green[pmo.lowoffdiag_begin[n] + idx] *= conj(ctem[index]);
                }
            }


            nstart += ct.block_dim[n];

        }

    }

}
