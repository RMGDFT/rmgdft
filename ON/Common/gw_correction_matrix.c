/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

   with GW correction, the eigenvalues are shifted, projected into localized orbital space.
   c

 */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"






void gw_correction_matrix(double *matS, double *Cij)
{
    int num_states = ct.num_states;
    int numst = ct.num_states;
    int ione = 1;    /* blas constants */

    double zero = 0., one = 1.;
    int st1;

    int i, j, li, lj, idx;
    /* for parallel libraries */
    int ictxt, mb, nb, npcol, nprow, llda, locc;
    int mycol, myrow;
    llda =  MXLLDA;
    locc =  MXLCOL;
    double *gw_shift, *correction_matrix00, *correction_matrix01;
    double xvec, xvec1;


    my_malloc(gw_shift, numst, double);
    my_malloc(correction_matrix00, llda * locc, double);
    my_malloc(correction_matrix01, llda * locc, double);


    for(st1 =0; st1 < numst/2; st1++)  gw_shift[st1] = 0.0;
    for(st1 =numst/2; st1 < numst; st1++)  gw_shift[st1] = 0.1;

    diag_eig_matrix(gamma_dis, gw_shift, pct.desca);



    /* work_dis = C * gamma_dis * C^*  */
    PSSYMM("r", "l", &numst, &numst, &one,
            gamma_dis, &ione, &ione, pct.desca,
            Cij, &ione, &ione, pct.desca, &zero, uu_dis, &ione, &ione, pct.desca);

    PSGEMM("N", "T", &numst, &numst, &numst, &one,
            uu_dis, &ione, &ione, pct.desca,
            Cij, &ione, &ione, pct.desca, &zero, work_dis, &ione, &ione, pct.desca);

    //work_dis2 = work_dis * S
    PSGEMM("N", "N", &numst, &numst, &numst, &one,
            work_dis, &ione, &ione, pct.desca,
            matS, &ione, &ione, pct.desca, &zero, work_dis2, &ione, &ione, pct.desca);
    // work_dis = S * work_dis1 , the final correction matrix 
    PSGEMM("N", "N", &numst, &numst, &numst, &one,
            matS, &ione, &ione, pct.desca,
            work_dis2, &ione, &ione, pct.desca, &zero, work_dis, &ione, &ione, pct.desca);

    // split work_dis, the correction matrix into 00 and 01, 
    // orbital overlap in the same unit cell Matrix 00  or by one unit cell translation in x
    // direction Matrix 01



    ictxt = pct.desca[1];
    mb = pct.desca[4];
    nb = pct.desca[5];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);



    double ax = get_celldm(0);
    for(li = 0; li < llda; li++)
    {
        for(lj = 0; lj < locc; lj++)
        {

            /*  li,lj are the  index of distributed matrix */
            /*  i,j are the  index of nondistributed matrix */
            i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
            j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

            idx = lj * llda + li;

            xvec = states[j].crds[0] - states[i].crds[0];
            xvec1 = states[i].crds[0]+ax - states[j].crds[0];
            if(xvec < ax/2.0 && xvec1 <ax/2.0) 
                error_handler(" ax is smaller than two-principle layers");
            else if(xvec < ax/2.0)
                correction_matrix00[idx] = work_dis[idx];
            else if(xvec < ax/2.0)
                correction_matrix01[idx] = work_dis[idx];
        }
    }

//  write out the matrix

	int fhand;

	int size, rank, ndims, gsizes[2], distribs[2];
	int order,  dargs[2], psizes[2];
	MPI_Datatype mpi_darray_lead;
	int mpi_amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
	MPI_File mpi_fhand ;
	MPI_Info fileinfo;
	MPI_Status status;
	MPI_Offset disp;



	size = nprow*npcol;
	ndims = 2;  
	distribs[0]= MPI_DISTRIBUTE_CYCLIC;
	distribs[1]= MPI_DISTRIBUTE_CYCLIC;
	dargs[0] = mb;
	dargs[1] = nb;
	psizes[0] = nprow;
	psizes[1] = npcol;
	order = MPI_ORDER_FORTRAN;

	gsizes[0] = ct.num_states;
	gsizes[1] = ct.num_states;
    rank = pct.gridpe;
	MPI_Type_create_darray(size, rank, ndims, gsizes, distribs, 
			dargs, psizes, order, MPI_DOUBLE, &mpi_darray_lead);
	MPI_Type_commit(&mpi_darray_lead);



    /* write out the matrices for left lead */

    MPI_Info_create(&fileinfo);
    char newname[MAX_PATH + 20];
	sprintf(newname, "%s%s", ct.outfile, ".gw_matrix");
    MPI_File_open(pct.grid_comm, newname, mpi_amode, fileinfo, &mpi_fhand);

    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, mpi_darray_lead, "native", MPI_INFO_NULL);
    idx = llda * locc;
    MPI_File_write_all(mpi_fhand, correction_matrix00, idx, MPI_DOUBLE, &status);
    MPI_File_write_all(mpi_fhand, correction_matrix01, idx, MPI_DOUBLE, &status);

    MPI_File_close(&mpi_fhand);
}


