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


#include "blas.h"
#include "transition.h"

template void DiagScalapack<double>(STATE *states, int numst, double *Hij_dist, double *Sij_dist);
template void DiagScalapack<std::complex<double>>(STATE *states, int numst, double *Hij_dist, double *Sij_dist);

template <typename KpointType>
void DiagScalapack(STATE *states, int numst, double *Hij_dist, double *Sij_dist)
{

    static Scalapack *MainSp = NULL;
    RmgTimer  *RT0 = new RmgTimer("3-DiagScalapack");
    if(MainSp == NULL)
    {
        int scalapack_groups = 1;
        int last = 1;
        MainSp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, numst,
                ct.scalapack_block_factor, last, pct.grid_comm);
    }

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    static std::vector<int> min_index;
    KpointType phase_k[27];
    std::complex<double> *phase_k_C = (std::complex<double> *)phase_k;

    if(min_index.size() == 0)
    {
        min_index.resize(MXLLDA * MXLCOL);
        MatrixKpointPhase (states, pct.desca, min_index);
    }


    bool participates = MainSp->Participates();

    int ione = 1;    /* blas constants */

    int info;
    double zero = 0., one = 1., alpha;
    int st1, st_g;
    double *eigs= new double[numst];

    std::vector<double> eigs_all, kweight, occ;
    eigs_all.resize(ct.num_kpts_pe * ct.num_states);
    occ.resize(ct.num_kpts_pe * ct.num_states);
    kweight.resize(ct.num_kpts_pe);
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        kweight[kpt] = ct.kp[pct.kstart + kpt].kweight;
        for(int st = 0; st < numst; st++)
        {
            occ[kpt * numst + st] = ct.kp[pct.kstart+kpt].kstate[st].occupation[0];
        }
    }

    /* for parallel libraries */
    int mb= pct.desca[4];
    int  mxllda2;
    mxllda2 = MXLLDA * MXLCOL;
    int mat_size = MXLLDA * MXLCOL * sizeof(double) * factor;

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

    KpointType *Hk = (KpointType *)uu_dis;
    KpointType *Sk = (KpointType *)l_s;
    KpointType *eigvec = (KpointType *)zz_dis;
    int *ipiv = new int[numst];


    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        if(ct.is_gamma)
        {
            memcpy(Hk, Hij_dist, mat_size);
            memcpy(eigvec, Hij_dist, mat_size);
            memcpy(Sk, Sij_dist, mat_size);
        }
        else
        {
            std::complex<double> I(0.0, 1.0);
            int idx = 0;
            double *kpt_xtal = ct.kp[pct.kstart + kpt].kpt;
            for(int ix = -1; ix <= 1; ix++)
            {
                for(int iy = -1; iy <= 1; iy++)
                {
                    for(int iz = -1; iz <= 1; iz++)
                    {
                        phase_k_C[idx] = std::exp(+I * (ix * kpt_xtal[0] + iy * kpt_xtal[1] + iz * kpt_xtal[2]));
                        idx++;
                    }
                }
            }
            for (int i = 0; i < MXLLDA * MXLCOL; i++)
            {
                Hk[i] = Hij_dist[i] * phase_k[ min_index[i] ];
                Sk[i] = Sij_dist[i] * phase_k[ min_index[i] ];
            }
            memcpy(eigvec, Hk, mat_size);
        }

        RmgTimer *RT1 = new RmgTimer("3-DiagScalapack: diag");
        if(participates)
            MainSp->generalized_eigenvectors(eigvec, Sk, eigs, Hk);

        MPI_Bcast(eigs, numst, MPI_DOUBLE, 0, pct.grid_comm);
        delete(RT1);

        for (st1 = 0; st1 < numst; st1++)
        {
            eigs_all[kpt * numst + st1 ] = eigs[st1];
        }




        RmgTimer *RT5 = new RmgTimer("3-DiagScalapack: pscal occ ");
        //   uu_dis = zz_dis *(occ_diag)
        memcpy(Hk, eigvec, mat_size);

        for(st1 = 0; st1 <  MXLCOL; st1++)
        {
                st_g = (st1/mb) * pct.scalapack_npcol * mb + pct.scalapack_mycol *mb +st1%mb;

                if(st_g >= numst) 
                    alpha = 0.0;
                else
                    alpha = occ[kpt * numst + st_g];

            for(int st2 = 0; st2 < MXLLDA; st2++)
            {
                Hk[st1 * MXLLDA + st2] *= alpha;
            }

        }

        delete(RT5);

        RmgTimer *RT3 = new RmgTimer("3-DiagScalapack: gemm ");

        if(ct.is_gamma)
        {
            pdgemm("N", "T", &numst, &numst, &numst, &one,
                    zz_dis, &ione, &ione, pct.desca,
                    uu_dis, &ione, &ione, pct.desca, &zero, mat_X, &ione, &ione, pct.desca);
        }
        else
        {
            std::complex<double> oneC(1.0), zeroC(0.0);
            pzgemm("N", "C", &numst, &numst, &numst, &oneC,
                    (std::complex<double> *)eigvec, &ione, &ione, pct.desca,
                    (std::complex<double> *)Hk, &ione, &ione, pct.desca, &zeroC, 
                    (std::complex<double> *)Sk, &ione, &ione, pct.desca);

            for(int idx = 0; idx <mxllda2; idx++) 
                mat_X[idx] = std::real(Sk[idx] * std::conj(phase_k[ min_index[idx] ]));

        }

        delete(RT3);


        RmgTimer *RT1b = new RmgTimer("3-DiagScalapack: (S^-1)H");

        /* Compute matrix theta = matB^-1 * Hij  */
        if(ct.is_gamma)
        {
            pdgetrf(&numst, &numst, (double *)Sij_dist, &ione, &ione, pct.desca, ipiv, &info);
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

        }
        else
        {
            for (int i = 0; i < MXLLDA * MXLCOL; i++)
            {
                Sk[i] = Sij_dist[i] * phase_k[ min_index[i] ];
                Hk[i] = Hij_dist[i] * phase_k[ min_index[i] ];
            }
            pzgetrf(&numst, &numst, (std::complex<double> *)Sk, &ione, &ione, pct.desca, ipiv, &info);
            if(info !=0)
            { 
                printf("\n error in pzgetrf in mg_eig.c INFO = %d\n", info);
                fflush(NULL);
                exit(0);
            }

            pzgetrs("N", &numst, &numst, (std::complex<double> *)Sk, &ione, &ione, pct.desca, ipiv, 
                    (std::complex<double> *)Hk, &ione, &ione, pct.desca, &info);
            double t1 = 2.0;
            //for(int i = 0; i < mxllda2; i++) uu_dis[i] = t1 * std::real(Hk[i] * phase_k[ min_index[i] ]);
            for(int i = 0; i < mxllda2; i++) uu_dis[i] = t1 * std::real(Hk[i]);
        }


        delete(RT1b);

        for(int st = 0; st < numst; st++)
        {
            ct.kp[pct.kstart+kpt].kstate[st].eig[0]= eigs_all[kpt * numst + st];
        }
    }

    if(pct.gridpe == 0) write_eigs(ct.kp[0].kstate);
    fflush(NULL);
    ct.efermi = Fill_on(eigs_all, kweight, occ, ct.occ_width, ct.nel, ct.occ_mix, ct.occ_flag, ct.mp_order);
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        for(int st = 0; st < numst; st++)
        {
            ct.kp[pct.kstart+kpt].kstate[st].occupation[0]= occ[kpt * numst + st];
        }
    }

    delete [] eigs;
    delete [] ipiv;
    delete(RT0);
}

