/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

   allocate memory for matrix in LCR part

L: left lead
C: conductor
R: right lead


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



void allocate_matrix_LCR ()
{
    int st, st1;
    long nbytes;
    char newname[MAX_PATH + 200];
    char msg[200];

    int idx, ntot, iprobe, idx_delta;
    int ixmin, ixmax;
    int pex, pey, pez;
    int NX, NY, NZ, PNX0, PNX1, PNY0, PNY1, PNZ0, PNZ1;
    int data_indicator;
    double x0_old, x_start, hx_new, hx_old, *potrho;
    int NX1, NX0, NX2, NY0, NY1, NY2, NZ0, NZ1, NZ2;
    int X0, X1, NX_start0, NX_start2;

    double tem;
    int ix, iy;


    int  ndim, i;


    /* Wait until everybody gets here */
    my_barrier ();

    ntot = 0;
    ndim = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
        ndim += ct.block_dim[i];
    }

    for (i = 1; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[i];
    }


    if (ndim != ct.num_states)
    {
        printf ("\n %d %d ndim not equaol to nC in read_data_part.c", ndim, ct.num_states);
        exit (0);
    }


    my_malloc_init( lcr[0].Stri, ntot, double );
    my_malloc_init( lcr[0].Htri, ntot, double );

    /* allocate memory for H00, H01, S00, S01  */
    for (iprobe =1; iprobe <= cei.num_probe; iprobe++)
    {
        idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].S00, idx, double );
        my_malloc_init( lcr[iprobe].S01, idx, double );
        my_malloc_init( lcr[iprobe].H00, idx, double );
        my_malloc_init( lcr[iprobe].H01, idx, double );

    }


    /* allocate memory for HLC, SLC  */
    for (iprobe =1; iprobe <= cei.num_probe; iprobe++)
    {
		i = cei.probe_in_block[iprobe - 1];
        idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].SCL, idx, double );
        my_malloc_init( lcr[iprobe].HCL, idx, double );

    }
	
	

    /*  allocate memory for density matrix  */

    my_malloc_init( lcr[0].density_matrix_tri, ntot, double );

    if (ct.runflag == 111 | ct.runflag == 112 | ct.runflag == 1121)
    {
        for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
        {
            lcr[iprobe].density_matrix_tri = lcr[0].density_matrix_tri;
        }
    }
    else
    {
        if(cei.probe_noneq == 0) 
        {
            for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
            {
                my_malloc_init( lcr[iprobe].density_matrix_tri, ntot, double );
                for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
                {	
                    my_malloc_init( lcr[iprobe].lcr_ne[idx_delta - 1].density_matrix_ne_tri, ntot, double );
                }
            }

        }

        else

        {

            my_malloc_init( lcr[1].density_matrix_tri, ntot, double );
            for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
            {	
                my_malloc_init( lcr[1].lcr_ne[idx_delta - 1].density_matrix_ne_tri, ntot, double );
            }

            for (iprobe = 2; iprobe <= cei.num_probe; iprobe++)
            {
                lcr[iprobe].density_matrix_tri = lcr[1].density_matrix_tri;
                for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
                {	
                    lcr[iprobe].lcr_ne[idx_delta - 1].density_matrix_ne_tri = 
                        lcr[1].lcr_ne[idx_delta - 1].density_matrix_ne_tri;
                }
            }
        }
    }



    int size  = (ct.state_end - ct.state_begin) * ct.num_states;
    if(size == 0 ) st = ct.num_states;
    my_malloc(Hij_00, size, double);
    my_malloc(Bij_00, size, double);



    fflush (NULL);


}    

