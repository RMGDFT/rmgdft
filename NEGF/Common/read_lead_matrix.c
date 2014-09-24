/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

   read_lead_matrix.c
   Functions to read matrix lcr[].H00, H01, S00, S01
   for lead only


 */



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void read_lead_matrix ()
{
    int iprobe, idx;
    int fhand;
    long nbytes, position;
    char *name;
    char newname[MAX_PATH + 200];
    int m,n, nprow, npcol, myrow, mycol, icrow, iccol;
    int i, j, k, ii, jj, ictxt, mb, nb, *desca;
    int nsize;
    double *temp, one = 1.0, zero = 0.0, *identity_matrix;
    int nL, ione = 1, dum;
    int idx_C, j1, jdiff;
    int t1, t2, jprobe;
    double *temp2;
    int idx1;
    double *tem_H00,  *tem_S00, *tem_H01, *tem_S01;
    double tem;

    int size, rank, ndims, gsizes[2], distribs[2];

    /* Wait until everyone gets here */
    my_barrier ();


    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;


    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    rank = myrow * pmo.ncol + mycol;

    size = pmo.nrow*pmo.ncol;
    ndims = 2;  


    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {


        name = lcr[iprobe].lead_name;
        sprintf (newname, "%s%s", name, ".matrix");


        idx = lcr[iprobe].num_states * lcr[iprobe].num_states;
        my_malloc_init( tem_H00, idx, double );
        my_malloc_init( tem_S00, idx, double );
        my_malloc_init( tem_H01, idx, double );
        my_malloc_init( tem_S01, idx, double );


        /* read the matrices for  leads */

        sprintf (newname, "%s%s", name, ".matrix");
        my_open( fhand, newname, O_RDWR, S_IREAD | S_IWRITE );
        nbytes = read(fhand, tem_H00, idx * sizeof(double));
        nbytes = read(fhand, tem_S00, idx * sizeof(double));
        nbytes = read(fhand, tem_H01, idx * sizeof(double));
        nbytes = read(fhand, tem_S01, idx * sizeof(double));
        if(nbytes != idx * sizeof(double) )
        {
            printf("\n end of file in lead matrix probe = %d\n", iprobe);
            exit(0);
        }
        close(fhand);

        /* For iprobe = 2, 4, 6, ... transpose H01, S01  */
        dum = 2 * (iprobe / 2);
        if(dum == iprobe) 
        {

            for(i = 0; i < lcr[iprobe].num_states; i++)
                for(j = i; j < lcr[iprobe].num_states; j++)
                {

                    idx  = i * lcr[iprobe].num_states + j;
                    idx1 = j * lcr[iprobe].num_states + i;

                    tem = tem_H01[idx]; 
                    tem_H01[idx] = tem_H01[idx1];
                    tem_H01[idx1] = tem;

                    tem = tem_S01[idx]; 
                    tem_S01[idx] = tem_S01[idx1];
                    tem_S01[idx1] = tem;
                }
        }

        if(iprobe == 3 && cei.num_probe == 3)
        {

            for(i = 0; i < lcr[iprobe].num_states; i++)
                for(j = i; j < lcr[iprobe].num_states; j++)
                {

                    idx  = i * lcr[iprobe].num_states + j;
                    idx1 = j * lcr[iprobe].num_states + i;

                    tem = tem_H01[idx]; 
                    tem_H01[idx] = tem_H01[idx1];
                    tem_H01[idx1] = tem;

                    tem = tem_S01[idx]; 
                    tem_S01[idx] = tem_S01[idx1];
                    tem_S01[idx1] = tem;
                }
        }


        /* now distribue the matrix into subset of processors */
        desca = &pmo.desc_lead[ (iprobe - 1) * DLEN];
        lead_mat_distribute(lcr[iprobe].H00, desca, tem_H00, iprobe);
        lead_mat_distribute(lcr[iprobe].S00, desca, tem_S00, iprobe);
        lead_mat_distribute(lcr[iprobe].H01, desca, tem_H01, iprobe);
        lead_mat_distribute(lcr[iprobe].S01, desca, tem_S01, iprobe);

        H01_to_HCL(tem_H01, lcr[iprobe].HCL, iprobe);
        H01_to_HCL(tem_S01, lcr[iprobe].SCL, iprobe);


        my_free(tem_H00);
        my_free(tem_S00);
        my_free(tem_H01);
        my_free(tem_S01);

    }

}
