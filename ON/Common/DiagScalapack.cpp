/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "make_conf.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"
#include "common_prototypes1.h"
#include "Scalapack.h"

//#include "main.h"
//#include "init_var.h"
//#include "RmgTimer.h"


#if !ELEMENTAL_LIBS
#include "blas.h"



void DiagScalapack(STATE *states, int numst, double *Hij_00, double *Bij_00, double *rho_matrix, double *theta_ptr)
{

    RmgTimer  *RT0 = new RmgTimer("3-DiagScalapack");
    int ione = 1;    /* blas constants */
    char *uplo = "u", *jobz = "v";

    int info;
    double zero = 0., one = 1.,half = 0.5,  alpha;
    int st1, st_g;
    double *eigs;

    /* for parallel libraries */
    int mb= pct.desca[4];
    int  mxllda2;
    mxllda2 = MXLLDA * MXLCOL;
    double *work1 = new double[mxllda2]();


    RmgTimer  *RT2 = new RmgTimer("3-DiagScalapack: cpdgemr2d");
    Cpdgemr2d(numst, numst, Hij_00, ione, ione, pct.descb, Hij, ione, ione,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Bij_00, ione, ione, pct.descb, matB, ione, ione,
            pct.desca, pct.desca[1]);

#if 0
dcopy(&mxllda2, matB, &ione, l_s, &ione);
PDTRAN(&numst, &numst, &half, l_s, &ione, &ione, pct.desca,
            &half, matB, &ione, &ione, pct.desca);

    
    pdlaset_( "A", &numst, &numst, &zero, &one, work1, &ione, &ione, pct.desca );
    int *ipiv = new int[2*numst]();
    pdgesv_(&numst, &numst, matB, &ione, &ione, pct.desca, ipiv, work1, &ione, &ione, pct.desca, &info);
    delete [] ipiv;

    pdgemm_("n", "n", &numst, &numst, &numst, &one,
                            work1, &ione, &ione, pct.desca, Hij, &ione,
                            &ione, pct.desca, &zero, uu_dis, &ione, &ione, pct.desca);
    dcopy(&mxllda2, uu_dis, &ione, Hij, &ione);
    PDTRAN(&numst, &numst, &half, Hij, &ione, &ione, pct.desca,
                &half, uu_dis, &ione, &ione, pct.desca);

    dcopy(&mxllda2, uu_dis, &ione, work1, &ione);
    delete(RT2);

Cpdgemr2d(numst, numst, Hij_00, ione, ione, pct.descb, Hij, ione, ione,
            pct.desca, pct.desca[1]);
Cpdgemr2d(numst, numst, Bij_00, ione, ione, pct.descb, matB, ione, ione,
            pct.desca, pct.desca[1]);
pdgemm_("n", "n", &numst, &numst, &numst, &one,
                        matB, &ione, &ione, pct.desca, work1, &ione,
                        &ione, pct.desca, &zero, uu_dis, &ione, &ione, pct.desca);
#endif

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
    dcopy(&mxllda2, Hij, &ione, uu_dis, &ione);
    dcopy(&mxllda2, matB, &ione, l_s, &ione);
//dcopy(&mxllda2, matB, &ione, l_s, &ione);
//PDTRAN(&numst, &numst, &half, l_s, &ione, &ione, pct.desca,
//           &half, matB, &ione, &ione, pct.desca);

    char *range = "a";
    double vx = 0.0;
    double tol = 0.0;
    int eigs_found, eigvs_found;
    double orfac = 0.0;
    int *iwork, *ifail, *iclustr, lwork;
    double *gap, lwork_tmp, *work2;
    int liwork_tmp, liwork;

    ifail = new int[numst];
    eigs = new double[numst];
    iclustr = new int[2 * pct.scalapack_nprow * pct.scalapack_npcol];
    gap = new double[pct.scalapack_nprow * pct.scalapack_npcol];

    lwork = -1;
    liwork = -1;

    pdsygvx (&ione, jobz, range, uplo, &numst, uu_dis, &ione, &ione, pct.desca,
            l_s, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
            &eigvs_found, eigs, &orfac, zz_dis, &ione, &ione, pct.desca, &lwork_tmp, &lwork,
            &liwork_tmp, &liwork, ifail, iclustr, gap, &info);

    if (info)
    {
        printf ("\n pdsygvx query failed, info is %d", info);
        exit(0);
    }


    /*set lwork and liwork */
    lwork = (int) lwork_tmp + 1;
    liwork = liwork_tmp;

    work2 = new double[lwork];
    iwork = new int[liwork];




    tol = 1.0e-15;
    pdsygvx (&ione, jobz, range, uplo, &numst, uu_dis, &ione, &ione, pct.desca,
            l_s, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
            &eigvs_found, eigs, &orfac, zz_dis, &ione, &ione, pct.desca, work2, &lwork,
            iwork, &liwork, ifail, iclustr, gap, &info);


    if (info)
    {
        printf ("\n pdsygvx failed, info is %d", info);
        exit(0);
    }



    delete(RT);


    RmgTimer *RT1 = new RmgTimer("3-DiagScalapack: calc_occ");
    for (st1 = 0; st1 < numst; st1++)
    {
        states[st1].eig[0] = eigs[st1];
    }

    delete [] ifail;
    delete [] iclustr;
    delete [] gap;
    delete [] work2;
    delete [] iwork;
    delete [] eigs;

    ct.efermi = fill_on(states, ct.occ_width, ct.nel, ct.occ_mix, numst, ct.occ_flag);
    delete(RT1);

    RmgTimer *RT5 = new RmgTimer("3-DiagScalapack: pscal occ ");
    //   uu_dis = zz_dis *(occ_diag)
    dcopy(&mxllda2, zz_dis, &ione, uu_dis, &ione);
    for(st1 = 0; st1 <  MXLCOL; st1++)
    {

        st_g = (st1/mb) * pct.scalapack_npcol * mb + pct.scalapack_mycol *mb +st1%mb;

        if(st_g >= numst) 
            alpha = 0.0;
        else
            alpha = states[st_g].occupation[0];
        dscal(&MXLLDA, &alpha, &uu_dis[st1 * MXLLDA], &ione);
    }

    delete(RT5);

    RmgTimer *RT3 = new RmgTimer("3-DiagScalapack: gemm ");

    pdgemm("N", "T", &numst, &numst, &numst, &one,
            uu_dis, &ione, &ione, pct.desca,
            zz_dis, &ione, &ione, pct.desca, &zero, mat_X, &ione, &ione, pct.desca);

    delete(RT3);


    RmgTimer *RT4 = new RmgTimer("3-DiagScalapack: XtoRow ");
    Cpdgemr2d(numst, numst, mat_X, ione, ione, pct.desca, rho_matrix, ione, ione,
            pct.descb, pct.desca[1]);
    delete(RT4);

    RmgTimer *RT1b = new RmgTimer("3-DiagScalapack: (S^-1)H");

    int *ipiv = new int[numst];
    /* Compute matrix theta = matB^-1 * Hij  */
    pdgetrf(&numst, &numst, matB, &ione, &ione, pct.desca, ipiv, &info);
    if(info !=0)
    { 
        printf("\n error in pdgetrf in mg_eig.c INFO = %d\n", info);
        fflush(NULL);
        exit(0);
    }

    dcopy(&mxllda2, Hij, &ione, uu_dis, &ione);
    pdgetrs("N", &numst, &numst, matB, &ione, &ione, pct.desca, ipiv, 
            uu_dis, &ione, &ione, pct.desca, &info);

    delete [] ipiv;


    double t1 = 2.0;
    dscal(&mxllda2, &t1, uu_dis, &ione);

    Cpdgemr2d(numst, numst, uu_dis, ione, ione, pct.desca, theta_ptr, ione, ione,
            pct.descb, pct.desca[1]);
    delete(RT1b);

    delete [] work1;


    delete(RT0);
}
#endif

