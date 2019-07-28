#include "negf_prototypes.h"
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
#include <complex>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"



void split_matrix_center ()
{

    int ntot, i, j, n, ii, jj, li, lj,nstart,idx, index;
    int *desca, ictxt, mb, nb, nprow, npcol, myrow, mycol;
    double distance[9], *Htem, *Stem, *Hyz[9], *Syz[9];
    double blength, clength, yvec, zvec;
    
    ntot = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
    }
    for (i = 1; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[i];
    }

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
        for(i = 0; i < 9; i++)
        {
            Hyz[i]   = &lcr[0].Htri_yz[pmo.diag_begin[n] + i * ntot];
            Syz[i]   = &lcr[0].Stri_yz[pmo.diag_begin[n] + i * ntot];
        }        


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

                Hyz[index][idx] = Htem[idx];
                Syz[index][idx] = Stem[idx];
            }
        }

        nstart += ct.block_dim[n];

    }



    nstart = 0;
    for(n = 0; n < ct.num_blocks -1; n++)
    {
        Htem  = &lcr[0].Htri[pmo.offdiag_begin[n]];
        Stem  = &lcr[0].Stri[pmo.offdiag_begin[n]];
        for(i = 0; i < 9; i++)
        {
            Hyz[i]   = &lcr[0].Htri_yz[pmo.offdiag_begin[n] + i * ntot];
            Syz[i]   = &lcr[0].Stri_yz[pmo.offdiag_begin[n] + i * ntot];
        }        


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

                Hyz[index][idx] = Htem[idx];
                Syz[index][idx] = Stem[idx];
            }
        }

        nstart += ct.block_dim[n];

    }


}
