/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "pmo.h"


void writeout_matrix_p ()
{

    int size, rank, ndims, gsizes[2], distribs[2];
    int order,  dargs[2], psizes[2];
    MPI_Datatype *mpi_darray_lead, *mpi_darray_diag, *mpi_darray_offdiag;
    int idx, i, iprobe;
    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;
    MPI_Info fileinfo;
    MPI_Status status;
    MPI_Offset disp;
    char newname[100];

    int ictxt = pmo.ictxt[pmo.myblacs];
    int mb = pmo.mblock;

    int nprow, npcol, myrow, mycol;

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    rank = myrow * pmo.ncol + mycol;

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

    if(rank == 0)
    {

        /* write out the matrices for leads */
        for(iprobe = 0; iprobe < cei.num_probe; iprobe++)
		{

			sprintf(newname, "%s%d%s", "lead_", iprobe, ".dat");
			MPI_Info_create(&fileinfo);
			MPI_File_open(COMM_EN2, newname, amode, fileinfo, &mpi_fhand);

			disp=0;
			MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, mpi_darray_lead[iprobe], "native", MPI_INFO_NULL);
			idx = pmo.mxllda_lead[iprobe]* pmo.mxlocc_lead[iprobe];
			MPI_File_write_all(mpi_fhand, lcr[iprobe + 1].H00, idx, MPI_DOUBLE, &status);
			MPI_File_write_all(mpi_fhand, lcr[iprobe + 1].S00, idx, MPI_DOUBLE, &status);
			MPI_File_write_all(mpi_fhand, lcr[iprobe + 1].H01, idx, MPI_DOUBLE, &status);
			MPI_File_write_all(mpi_fhand, lcr[iprobe + 1].S01, idx, MPI_DOUBLE, &status);

			MPI_File_close(&mpi_fhand);
        }
		

        /* write the H00 and S00 for center part, in the form of tri-diagonal blocks */
        /* the file names will be center_diag_#.dat and center_offdiag_#.dat 
         * */

        for(i =0; i < ct.num_blocks; i++)
        {

            sprintf(newname, "%s%d%s", "center_diag_", i, ".dat");

            MPI_Info_create(&fileinfo);
            MPI_File_open(COMM_EN2, newname, amode, fileinfo, &mpi_fhand);

            disp=0;
            MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, mpi_darray_diag[i], "native", MPI_INFO_NULL);
            idx = pmo.mxllda_cond[i]* pmo.mxlocc_cond[i];
            MPI_File_write_all(mpi_fhand, &lcr[0].Htri[pmo.diag_begin[i]], idx, MPI_DOUBLE, &status);
            MPI_File_write_all(mpi_fhand, &lcr[0].Stri[pmo.diag_begin[i]], idx, MPI_DOUBLE, &status);
            MPI_File_close(&mpi_fhand);
        }

        for(i =0; i < ct.num_blocks-1; i++)
        {
            sprintf(newname, "%s%d%s", "center_offdiag_", i, ".dat");
            MPI_Info_create(&fileinfo);
            MPI_File_open(COMM_EN2, newname, amode, fileinfo, &mpi_fhand);

            disp=0;
            MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, mpi_darray_offdiag[i], "native", MPI_INFO_NULL);
            idx = pmo.mxllda_cond[i]* pmo.mxlocc_cond[i+1];
            MPI_File_write_all(mpi_fhand, &lcr[0].Htri[pmo.offdiag_begin[i]], idx, MPI_DOUBLE, &status);
            MPI_File_write_all(mpi_fhand, &lcr[0].Stri[pmo.offdiag_begin[i]], idx, MPI_DOUBLE, &status);
            MPI_File_close(&mpi_fhand);

        }



    }

    if (pct.gridpe == 0)
        printf ("\n THE MATRICES are written out ");
}


