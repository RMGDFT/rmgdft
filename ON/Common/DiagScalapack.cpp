/************************** SVN Revision Information **************************
 **    $Id$    **
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
#include "prototypes_on.h"

//#include "main.h"
//#include "init_var.h"
//#include "RmgTimer.h"


#if !ELEMENTAL_LIBS
#include "blas.h"

Scalapack *MainSp;
extern Scalapack *MainSp;



void DiagScalapack(STATE *states, int numst, double *Hij_dist, double *Sij_dist)
{

    RmgTimer  *RT0 = new RmgTimer("3-DiagScalapack");
    if(!MainSp) 
    {
        int scalapack_groups = 1;
        int last = 1;
        MainSp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, numst,
                ct.scalapack_block_factor, last, pct.grid_comm);
    }

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    bool participates = MainSp->Participates();

    int ione = 1;    /* blas constants */

    int info;
    double zero = 0., one = 1., alpha;
    int st1, st_g;
    double *eigs= new double[numst];

    /* for parallel libraries */
    int mb= pct.desca[4];
    int  mxllda2;
    mxllda2 = MXLLDA * MXLCOL;
    int mat_size = MXLLDA * MXLCOL * sizeof(double) * factor;

    RmgTimer *RT = new RmgTimer("3-DiagScalapack: pdsygvx ");
    /* If I'm in the process grid, execute the program */
    if (pct.scalapack_myrow < 0)
    {  
        printf("\n nprow, npcol %d %d", pct.scalapack_nprow, pct.scalapack_npcol);
        printf("\n we should use all proc for diag. somthing wrong");
        printf("\n gridpe = %d", pct.gridpe);
        exit(0);
    }


    /* 
     * SOLVE THE GENERALIZED EIGENVALUE PROBLEM:  m * z = lambda * matS * z 
     */

    /* Transform the generalized eigenvalue problem to a sStandard form */
    memcpy(uu_dis, Hij_dist, mat_size);
    memcpy(zz_dis, Hij_dist, mat_size);
    memcpy(l_s, Sij_dist, mat_size);

    if(participates)
        MainSp->generalized_eigenvectors(zz_dis, l_s, eigs, uu_dis);

    delete(RT);


    MPI_Bcast(eigs, numst, MPI_DOUBLE, 0, pct.grid_comm);

    RmgTimer *RT1 = new RmgTimer("3-DiagScalapack: calc_occ");
    for (st1 = 0; st1 < numst; st1++)
    {
        states[st1].eig[0] = eigs[st1];
    }

    delete [] eigs;
    if(ct.spin_flag)
    {
        get_opposite_eigvals( states );
    }
    if(pct.gridpe == 0) write_eigs(states);
    fflush(NULL);
    /* Generate new density */

    ct.efermi = Fill_on(states, ct.occ_width, ct.nel, ct.occ_mix, numst, ct.occ_flag, ct.mp_order);


    delete(RT1);

    RmgTimer *RT5 = new RmgTimer("3-DiagScalapack: pscal occ ");
    //   uu_dis = zz_dis *(occ_diag)
    memcpy(uu_dis, zz_dis, mat_size);
    
    for(st1 = 0; st1 <  MXLCOL; st1++)
    {

        st_g = (st1/mb) * pct.scalapack_npcol * mb + pct.scalapack_mycol *mb +st1%mb;

        if(st_g >= numst) 
            alpha = 0.0;
        else
            alpha = states[st_g].occupation[0];

        for(int st2 = 0; st2 < MXLLDA; st2++)
            uu_dis[st1 * MXLLDA + st2] *= alpha;
                 
    }

    delete(RT5);

    RmgTimer *RT3 = new RmgTimer("3-DiagScalapack: gemm ");

    pdgemm("N", "T", &numst, &numst, &numst, &one,
            uu_dis, &ione, &ione, pct.desca,
            zz_dis, &ione, &ione, pct.desca, &zero, mat_X, &ione, &ione, pct.desca);

    delete(RT3);


    RmgTimer *RT1b = new RmgTimer("3-DiagScalapack: (S^-1)H");

    int *ipiv = new int[numst];
    /* Compute matrix theta = matB^-1 * Hij  */
    pdgetrf(&numst, &numst, Sij_dist, &ione, &ione, pct.desca, ipiv, &info);
    if(info !=0)
    { 
        printf("\n error in pdgetrf in mg_eig.c INFO = %d\n", info);
        fflush(NULL);
        exit(0);
    }

    memcpy(uu_dis, Hij_dist, mat_size);
    pdgetrs("N", &numst, &numst, Sij_dist, &ione, &ione, pct.desca, ipiv, 
            uu_dis, &ione, &ione, pct.desca, &info);

    double t1 = 2.0;
    for(int i = 0; i < mxllda2; i++) uu_dis[i] *= t1;

    int lwork = -1;
    int liwork = -1;
    double lwork_tmp;
    int liwork_tmp;
    pdgetri(&numst, Sij_dist, &ione, &ione, pct.desca, ipiv, &lwork_tmp, &lwork, &liwork_tmp, &liwork, &info);
    lwork = (int)lwork_tmp + 8;
    liwork = (int)liwork_tmp + 8;
    double *nwork = new double[lwork];
    int *iwork = new int[liwork];
    pdgetri(&numst, Sij_dist, &ione, &ione, pct.desca, ipiv, nwork, &lwork, iwork, &liwork, &info);


    delete(RT1b);

    delete [] nwork;
    delete [] iwork;
    delete [] ipiv;


    delete(RT0);
}
#endif

