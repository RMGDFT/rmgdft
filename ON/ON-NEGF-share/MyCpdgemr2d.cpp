
/************************** SVN Revision Information **************************
 **    $Id: DiagScalapack.cpp 4470 2018-08-07 20:48:04Z luw $    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"
#include "common_prototypes1.h"
#include "Scalapack.h"
#include "blacs.h"

void mat_global_to_dist(double *a_dist, int *desca, double *a_glob);
void mat_dist_to_global(double *a_dist, int *desca, double *a_glob);
void MyCpdgemr2d(int M,int N, double *A, int *desca, double *B, int *descb)
{

    if(ct.use_cpdgemr2d)
    {
        int ione = 1;
        Cpdgemr2d(M, N, A, ione, ione, desca, B, ione, ione, descb, desca[1]);
    }
    else
    {
        size_t msize = M * N;
        double *A_glob = new double[msize];

        mat_dist_to_global(A, desca, A_glob);

        mat_global_to_dist(B, descb, A_glob);

        delete [] A_glob;
    }
}

void mat_global_to_dist(double *a_dist, int *desca, double *a_glob)
{
    int mycol, myrow, nprow, npcol;
    int ictxt=desca[1], m = desca[2], n = desca[3], mb=desca[4], nb=desca[5], mxllda = desca[8];

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int izero = 0;
    int mxli = numroc(&m, &mb, &myrow, &izero, &nprow);
    int mxlocc = numroc(&n, &nb, &mycol, &izero, &npcol);


    for(int j =0; j < mxli; j++)
    {
        for(int k=0; k < mxlocc; k++)
        {

            /*  distributed local index (j,k) is (jj, kk) in global matrix
             */

            int jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb;
            int kk = (k/nb ) * npcol * nb + mycol * nb + k - k/nb * nb;

            a_dist[k * mxllda + j] = a_glob[kk * m  + jj] ;

        }
    }

}


void mat_dist_to_global(double *a_dist, int *desca, double *a_glob)
{
    int mycol, myrow, nprow, npcol;
    int ictxt=desca[1], m = desca[2], n = desca[3], mb=desca[4], nb=desca[5], mxllda = desca[8];

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int izero = 0;
    int mxli = numroc(&m, &mb, &myrow, &izero, &nprow);
    int mxlocc = numroc(&n, &nb, &mycol, &izero, &npcol);


    for(int i = 0; i < m *n; i++) a_glob[i] = 0.0;

    for(int j =0; j < mxli; j++)
    {
        for(int k=0; k < mxlocc; k++)
        {

            /*  distributed local index (j,k) is (jj, kk) in global matrix
             */

            int jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb;
            int kk = (k/nb ) * npcol * nb + mycol * nb + k - k/nb * nb;

            a_glob[kk * m  + jj] = a_dist[k * mxllda + j];

        }
    }

    int length = m*n;
    MPI_Allreduce(MPI_IN_PLACE, a_glob, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
}


