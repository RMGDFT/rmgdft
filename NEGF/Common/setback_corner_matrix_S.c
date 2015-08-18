/************************** SVN Revision Information **************************
 **    $Id: 
******************************************************************************/
 


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


/* 
* The matrix elements in the corners  (the related orbitals are in the
* boarder of the leads ) need to be same as in the leads. In our matrix
* calculations, peoridcal condition is used even when the left and right
* leads, up and down lead are totally different. The related matrix
* element is artificially changed. This function copy back from lead
* matrix. 
* 
 */


void  setback_corner_matrix_S()
{
    int i, j, ii, jj, iii, jjj, li, lj, maxli;
    int iistart, jjstart, limb, ljnb;
    int mycol, myrow, nprow, npcol;
    int ictxt, mb, nb, *desca, n0;
    int idx_local, idx_global, item;
    int iprobe, jprobe, idx, ni,  nL, nR;
    int ni0, ni1;
    int nmax,idx1, nlead0, n2;

    double *temp;

    Cblacs_gridinfo (pmo.ictxt[pmo.myblacs], &nprow, &npcol, &myrow, &mycol);


    mb = pmo.mblock;
    nb = pmo.mblock;

    nmax = 0;
    for(iprobe = 1; iprobe <= cei.num_probe ; iprobe++)
    {

        nmax = rmg_max (nmax, lcr[iprobe].num_states);
    }
    
    my_malloc(temp, nmax * nmax, double);


    for(iprobe = 1; iprobe <= cei.num_probe ; iprobe++)
    {

        /* copy distributed S00 matrix to global temp */

        for(idx = 0; idx < nmax * nmax; idx++) temp[idx] = 0.0;

        for (lj = 0; lj < pmo.mxlocc_lead[iprobe-1]; lj++)
        {
            for (li = 0; li < pmo.mxllda_lead[iprobe-1]; li++)
            {
                i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;
                
                idx = j * lcr[iprobe].num_states + i;
                idx1 = lj * pmo.mxllda_lead[iprobe-1] + li;
                
                temp[idx] = lcr[iprobe].S00[idx1];
            }
        }
        
        n2 = nmax * nmax;
        comm_sums(temp, &n2, COMM_EN2);
        /* now temp = S00, not distributed */

        ni = cei.probe_in_block[iprobe-1];  /* n0 => block index */

        /* determine which part of matrix should be modified */

        if(iprobe%2 == 0 ) 
        { 
            ni0 = ct.block_dim[ni] - lcr[iprobe].num_states/2;
            ni1 = ct.block_dim[ni];
            nlead0 = ct.block_dim[ni] - lcr[iprobe].num_states;
        }
        else
        {
            ni0 = 0;
            ni1 = lcr[iprobe].num_states/2;
            nlead0 = 0;
        }

        if(iprobe == 3 && cei.num_probe == 3)
        {
            ni0 = lcr[iprobe].num_states/2;
            ni1 =  lcr[iprobe].num_states;
            nlead0 = 0;
        }



        for (lj = 0; lj < pmo.mxlocc_cond[ni]; lj++)
        {
            for (li = 0; li < pmo.mxllda_cond[ni]; li++)
            {
                i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;


                /* set the lead image interaction to be zero */

                idx = lj * pmo.mxllda_cond[ni] + li;
                idx1 = (j- nlead0) * lcr[iprobe].num_states + i - nlead0;

                if( i >= ni0 && i < ni1 && j >=ni0 && j < ni1)
                {
                    lcr[0].Stri[pmo.diag_begin[ni] + idx] = temp[idx1];
                }


            }  /* end li */
        } /* end lj */

    } 

    my_free(temp);

}
