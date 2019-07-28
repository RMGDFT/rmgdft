#include "negf_prototypes.h"
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


/* set the image interaction of the leads to be zero 
*  when two leads are connected to the same block or neighbor blocks, 
*  their image interaction must set to be zero due to we calculate 
* the matrix with image interaction included.
 */


void  zero_lead_image(double *tri)
{
    int i, j, ii, jj, iii, jjj, li, lj, maxli;
    int iistart, jjstart, limb, ljnb;
    int mycol, myrow, nprow, npcol;
    int ictxt, mb, nb, *desca, n0;
    int idx_local, idx_global, item;
    int iprobe, jprobe, idx, ni, nj, nL, nR;
    int ni0, ni1, nj0, nj1;
    Cblacs_gridinfo (pmo.ictxt[pmo.myblacs], &nprow, &npcol, &myrow, &mycol);

    mb = pmo.mblock;
    nb = pmo.mblock;

    for(iprobe = 1; iprobe <= cei.num_probe ; iprobe++)
    {
        for(jprobe = iprobe+1; jprobe <= cei.num_probe ; jprobe++)
        {


            ni = cei.probe_in_block[iprobe-1];  /* n0 => block index */
            nj = cei.probe_in_block[jprobe-1];  /* n0 => block index */

            if(ni == nj)  /* two leads connect to the same block */
            {

                if(iprobe%2 == 1) 
                {

                    nL = lcr[iprobe].num_states;
                    nR = lcr[jprobe].num_states; 
                }
                else if( jprobe%2 == 1)
                {
                    nR = lcr[iprobe].num_states;
                    nL = lcr[jprobe].num_states; 
                }
                else
                {
                    printf("\n two probes at the same block must be one even and one odd \n");
                    printf("\n iprobe = %d, jprobe = %d\n", iprobe, jprobe);

                    exit(0);
                }


                /* li, lj: the index in distributed matrix
                 * i, j: related index in the global matrix
                 */


                for (lj = 0; lj < pmo.mxlocc_cond[ni]; lj++)
                {
                    for (li = 0; li < pmo.mxllda_cond[ni]; li++)
                    {
                        i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                        j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;


                        /* set the lead image interaction to be zero */

                        idx = lj * pmo.mxllda_cond[ni] + li;
                        if( j < nL && i >= ct.block_dim[ni] - nR ) 
                        {
                            tri[pmo.diag_begin[ni] + idx] = 0;
                        }

                        if( i < nL && j >= ct.block_dim[ni] - nR )  
                        {

                            tri[pmo.diag_begin[ni] + idx] = 0;
                        }

                    }  /* end li */
                } /* end lj */

            } /* end if (ni == nj) */


            /* determine where the image  is */
            if(ni%2 == 0 ) 
            { 
                ni0 = ct.block_dim[ni] - lcr[iprobe].num_states;
                ni1 = ct.block_dim[ni];
            }
            else
            {
                ni0 = 0;
                ni1 = lcr[iprobe].num_states;
            }
            if(nj%2 == 0 ) 
            { 
                nj0 = ct.block_dim[nj] - lcr[jprobe].num_states;
                nj1 = ct.block_dim[nj];
            }
            else
            {
                nj0 = 0;
                nj1 = lcr[jprobe].num_states;
            }

            if( nj - ni == 1)  /* two leads are connected to the neighboring blocks */
            {

                for (lj = 0; lj < pmo.mxlocc_cond[nj]; lj++)
                {
                    for (li = 0; li < pmo.mxllda_cond[ni]; li++)
                    {
                        i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                        j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;


                        /* set the lead image interaction to be zero */

                        idx = lj * pmo.mxllda_cond[ni] + li;
                        if( i >= ni0 && i < ni1 && j >=nj0 && j < nj1)
                        {
                            tri[pmo.offdiag_begin[ni] + idx] = 0;
                        }


                    }  /* end li */
                } /* end lj */
                
            } /* end if ( nj - ni == 1) */

            if( ni - nj == 1)  /* two leads are connected to the neighboring blocks */
            {

                for (lj = 0; lj < pmo.mxlocc_cond[ni]; lj++)
                {
                    for (li = 0; li < pmo.mxllda_cond[nj]; li++)
                    {
                        i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
                        j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;


                        /* set the lead image interaction to be zero */

                        idx = lj * pmo.mxllda_cond[nj] + li;
                        if( i >= nj0 && i < nj1 && j >=ni0 && j < ni1)
                        {
                            tri[pmo.offdiag_begin[nj] + idx] = 0;
                        }


                    }  /* end li */
                } /* end lj */
                
            } /* end if ( ni - nj == 1) */

        }

    }
}

