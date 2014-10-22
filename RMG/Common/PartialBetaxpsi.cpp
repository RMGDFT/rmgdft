/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"


template void PartialBetaxpsi<double> (int ion, fftw_plan p2, double * newsintR_x, double * newsintR_y,
                       double * newsintR_z, double * newsintI_x, double * newsintI_y, double * newsintI_z,
                       ION * iptr, Kpoint<double> **Kptr);
template void PartialBetaxpsi<std::complex<double>> (int ion, fftw_plan p2, double * newsintR_x, double * newsintR_y,
                       double * newsintR_z, double * newsintI_x, double * newsintI_y, double * newsintI_z,
                       ION * iptr,Kpoint<std::complex<double>> **Kptr);

template <typename OrbitalType> void PartialBetaxpsi (int ion, fftw_plan p2, double * newsintR_x, double * newsintR_y,
                       double * newsintR_z, double * newsintI_x, double * newsintI_y, double * newsintI_z,
                       ION * iptr,Kpoint<OrbitalType> **Kptr)
{

    int idx, kidx, istate, size, nh, index;
    int prjcount, count;
    int incx = 1, *pidx, ip;
    double *workR;
    double *beta_x, *beta_y, *beta_z;
    double *temp_psiR;
    SPECIES *sp;
    int P0_BASIS;
    double *workI, *pI, *temp_psiI, *pR;

    nh = ct.sp[ct.ions[ion].species].nh;
    sp = &ct.sp[ct.ions[ion].species];
    P0_BASIS = get_P0_BASIS();

    count = pct.idxptrlen[ion];
    pidx = pct.nlindex[ion];

    if (!count) return;

    workR = new double[2* P0_BASIS];
    workI = workR + P0_BASIS;


    size = nh * P0_BASIS;

    beta_x = new double[ 3 * size]; 

    beta_y = beta_x + size;
    beta_z = beta_y + size;

    for (idx = 0; idx < size; idx++)
    {
        beta_x[idx] = ZERO;
        beta_y[idx] = ZERO;
        beta_z[idx] = ZERO;
    }

    /*partial_beta(ion, beta_x, beta_y,beta_z, iptr, p1, p2); */
    get_derweight (ion, beta_x, beta_y, beta_z, iptr, p2);


    for (kidx = 0; kidx < ct.num_kpts; kidx++)
    {


        if(!ct.is_gamma)
        {
            pR = pct.phaseptr[ion];
            pR += 2 * kidx * sp->nldim * sp->nldim * sp->nldim;
            pI = pR + sp->nldim * sp->nldim * sp->nldim;
        }

        for (istate = 0; istate < ct.num_states; istate++)
        {


            /*Gather_psi is not necessary, getting pointer should be enough
              gather_psi(temp_psiR, temp_psiI, sta, 0); */

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                if(ct.is_gamma)
                    workR[idx] = std::real(Kptr[kidx]->Kstates[istate].psi[idx]);
                else
                {
                    workR[idx] = std::real(Kptr[kidx]->Kstates[istate].psi[idx]) * pR[idx] - std::imag(Kptr[kidx]->Kstates[istate].psi[idx]) * pI[idx];
                    workI[idx] = std::imag(Kptr[kidx]->Kstates[istate].psi[idx]) * pR[idx] + std::real(Kptr[kidx]->Kstates[istate].psi[idx]) * pI[idx];
                }
            }

            index =
                ion * ct.num_kpts * ct.num_states * ct.max_nl +
                kidx * ct.num_states * ct.max_nl + istate * ct.max_nl;
            for (ip = 0; ip < nh; ip++)
            {

                newsintR_x[index+ip] =
                    get_vel() * QMD_ddot (P0_BASIS, workR, incx, &beta_x[ip*P0_BASIS], incx);
                newsintR_y[index+ip] =
                    get_vel() * QMD_ddot (P0_BASIS, workR, incx, &beta_y[ip*P0_BASIS], incx);
                newsintR_z[index+ip] =
                    get_vel() * QMD_ddot (P0_BASIS, workR, incx, &beta_z[ip*P0_BASIS], incx);
                if(!ct.is_gamma)
                {
                    newsintI_x[index+ip] =
                        get_vel() * QMD_ddot (P0_BASIS, workI, incx, &beta_x[ip*P0_BASIS], incx);
                    newsintI_y[index+ip] =
                        get_vel() * QMD_ddot (P0_BASIS, workI, incx, &beta_y[ip*P0_BASIS], incx);
                    newsintI_z[index+ip] =
                        get_vel() * QMD_ddot (P0_BASIS, workI, incx, &beta_z[ip*P0_BASIS], incx);
                }

            }               /*end for ip */
        }                   /*end if count */
    }                       /* end for istate */

    delete[] workR;
    delete[] beta_x;
}
