#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"



void read_matrix_pp ()
{

    int size, rank, ndims, gsizes[2], distribs[2];
    int order,  dargs[2], psizes[2];
    MPI_Datatype *mpi_darray_lead, *mpi_darray_diag, *mpi_darray_offdiag;
    int idx, i, iprobe;
    int amode = MPI_MODE_RDWR;
    MPI_File mpi_fhand ;
    MPI_Info fileinfo;
    MPI_Status status;
    MPI_Offset disp;
    char newname[100];

    int ictxt = pmo.ictxt[pmo.myblacs];
    int fhand, *desca;
    size_t nbytes;
    double *tem_H00, *tem_S00, *tem_H01, *tem_S01;

    int nprow, npcol, myrow, mycol;

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    my_malloc(mpi_darray_lead, cei.num_probe+1, MPI_Datatype);
    my_malloc(mpi_darray_diag, ct.num_blocks, MPI_Datatype);
    my_malloc(mpi_darray_offdiag, ct.num_blocks, MPI_Datatype);
    size = pmo.nrow*pmo.ncol;
    ndims = 2;  
    distribs[0]= MPI_DISTRIBUTE_CYCLIC;
    distribs[1]= MPI_DISTRIBUTE_CYCLIC;
    dargs[0] = pmo.mblock;
    dargs[1] = pmo.mblock;
    psizes[0] = pmo.nrow;
    psizes[1] = pmo.ncol;
    order = MPI_ORDER_FORTRAN;

    rank = myrow *pmo.ncol + mycol ;
    
    for(iprobe = 0; iprobe < cei.num_probe; iprobe++)
    {
        gsizes[0] = lcr[iprobe+1].num_states;
        gsizes[1] = lcr[iprobe+1].num_states;

        MPI_Type_create_darray(size, rank, ndims, gsizes, distribs, 
                dargs, psizes, order, MPI_DOUBLE, &mpi_darray_lead[iprobe]);
        MPI_Type_commit(&mpi_darray_lead[iprobe]);
    }

    for(i = 0; i < ct.num_blocks; i++)
    {
        gsizes[0] = ct.block_dim[i];
        gsizes[1] = ct.block_dim[i];
        MPI_Type_create_darray(size, rank, ndims, gsizes, distribs, 
                dargs, psizes, order, MPI_DOUBLE, &mpi_darray_diag[i]);
        MPI_Type_commit(&mpi_darray_diag[i]);
    }

    for(i = 0; i < ct.num_blocks-1; i++)
    {
        gsizes[0] = ct.block_dim[i];
        gsizes[1] = ct.block_dim[i+1];
        MPI_Type_create_darray(size, rank, ndims, gsizes, distribs, 
                dargs, psizes, order, MPI_DOUBLE, &mpi_darray_offdiag[i]);
        MPI_Type_commit(&mpi_darray_offdiag[i]);
    }

      MPI_Comm_rank(COMM_EN1, &rank);




    /* read the matrices for leads */

    for(iprobe = 0; iprobe < cei.num_probe; iprobe++)
	{	

        idx = lcr[iprobe + 1].num_states * lcr[iprobe + 1].num_states;
        my_malloc_init( tem_H00, idx, double );
        my_malloc_init( tem_S00, idx, double );
        my_malloc_init( tem_H01, idx, double );
        my_malloc_init( tem_S01, idx, double );
                                                                                                                                 
        /* read the matrices for  leads */

		sprintf(newname, "%s%s%d%s", pct.image_path[pct.thisimg], "lead_", iprobe, ".dat");
        fhand = open(newname, O_RDWR);
        nbytes = read(fhand, tem_H00, idx * sizeof(double));
        nbytes = read(fhand, tem_S00, idx * sizeof(double));
        nbytes = read(fhand, tem_H01, idx * sizeof(double));
        nbytes = read(fhand, tem_S01, idx * sizeof(double));
        if(nbytes != (size_t)idx * sizeof(double) )
        {
            printf("\n end of file in lead matrix probe = %d\n", iprobe + 1);
            exit(0);
        }
        close(fhand);


        /* now distribue the matrix into subset of processors */
        desca = &pmo.desc_lead[ (iprobe) * DLEN];
        lead_mat_distribute(lcr[iprobe + 1].H00, desca, tem_H00, iprobe + 1);
        lead_mat_distribute(lcr[iprobe + 1].S00, desca, tem_S00, iprobe + 1);
        lead_mat_distribute(lcr[iprobe + 1].H01, desca, tem_H01, iprobe + 1);
        lead_mat_distribute(lcr[iprobe + 1].S01, desca, tem_S01, iprobe + 1);


        H01_to_HCL(tem_H01, lcr[iprobe + 1].HCL, iprobe + 1);
        H01_to_HCL(tem_S01, lcr[iprobe + 1].SCL, iprobe + 1);


        my_free(tem_H00);
        my_free(tem_S00);
        my_free(tem_H01);
        my_free(tem_S01);

		if(ct.runflag == 100) return;
    }

/* ------------------------------------------------------------ */


	/* readthe H00 and S00 for center part, in the form of tri-diagonal blocks */

    for(i =0; i < ct.num_blocks; i++)
    {

        sprintf(newname, "%s%s%d%s", pct.image_path[pct.thisimg], "center_diag_", i, ".dat");

        MPI_Info_create(&fileinfo);
        MPI_File_open(COMM_EN2, newname, amode, fileinfo, &mpi_fhand);

        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, mpi_darray_diag[i], "native", MPI_INFO_NULL);
        idx = pmo.mxllda_cond[i]* pmo.mxlocc_cond[i];
        MPI_File_read(mpi_fhand, &lcr[0].Htri[pmo.diag_begin[i]], idx, MPI_DOUBLE, &status);
        MPI_File_read(mpi_fhand, &lcr[0].Stri[pmo.diag_begin[i]], idx, MPI_DOUBLE, &status);
        MPI_File_close(&mpi_fhand);

    }

    for(i =0; i < ct.num_blocks-1; i++)
    {
        sprintf(newname, "%s%s%d%s", pct.image_path[pct.thisimg], "center_offdiag_", i, ".dat");

        MPI_Info_create(&fileinfo);
        MPI_File_open(COMM_EN2, newname, amode, fileinfo, &mpi_fhand);

        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, mpi_darray_offdiag[i], "native", MPI_INFO_NULL);
        idx = pmo.mxllda_cond[i]* pmo.mxlocc_cond[i+1];
        MPI_File_read(mpi_fhand, &lcr[0].Htri[pmo.offdiag_begin[i]], idx, MPI_DOUBLE, &status);
        MPI_File_read(mpi_fhand, &lcr[0].Stri[pmo.offdiag_begin[i]], idx, MPI_DOUBLE, &status);
        MPI_File_close(&mpi_fhand);

    }



    if (pct.gridpe == 0)
        printf ("\n THE MATRICES are read in ");
}
