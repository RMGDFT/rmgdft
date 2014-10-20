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
    int alloc, prjcount, count;
    int incx = 1, *pidx, ip;
    double *workR;
    double *beta_x, *beta_y, *beta_z;
    double *temp_psiR;
    STATE *sta;
#if !GAMMA_PT
    double *workI, *pI, *temp_psiI, *pR;
#endif

    nh = ct.sp[ct.ions[ion].species].nh;
    alloc = get_P0_BASIS();
    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;

    count = pct.idxptrlen[ion];
    pidx = pct.nlindex[ion];

    if (count)
    {
#if GAMMA_PT
        workR = new double[ alloc];
#else
        workR = new double[2* alloc];
        workI = workR + alloc;
#endif
    }


    size = nh * pct.idxptrlen[ion];
    if (size)
    {
#if !FDIFF_BETA
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
#else

        beta_x = pct.weight_derx[ion];
        beta_y = pct.weight_dery[ion];
        beta_z = pct.weight_derz[ion];

#endif
    }
    else
    {
        beta_x = NULL;
        beta_y = NULL;
        beta_z = NULL;
    }

    for (kidx = 0; kidx < ct.num_kpts; kidx++)
    {


#if !GAMMA_PT
        pR = pct.phaseptr[ion];
        pR += 2 * kidx * count;
        pI = pR + count;
#endif

        for (istate = 0; istate < ct.num_states; istate++)
        {

            if (count)
            {

                /*Gather_psi is not necessary, getting pointer should be enough
                   gather_psi(temp_psiR, temp_psiI, sta, 0); */

                for (idx = 0; idx < count; idx++)
                {
#if GAMMA_PT
                    workR[idx] = std::real(Kptr[kidx]->Kstates->psi[pidx[idx]]);
#else
                    workR[idx] = std::real(Kptr[kidx]->Kstates->psi[pidx[idx]]) * pR[idx] - std::imag(Kptr[kidx]->Kstates->psi[pidx[idx]]) * pI[idx];
                    workI[idx] = std::imag(Kptr[kidx]->Kstates->psi[pidx[idx]]) * pR[idx] + std::real(Kptr[kidx]->Kstates->psi[pidx[idx]]) * pI[idx];
#endif
                }

                prjcount = 0;
                index =
                    ion * ct.num_kpts * ct.num_states * ct.max_nl +
                    kidx * ct.num_states * ct.max_nl + istate * ct.max_nl;
                for (ip = 0; ip < nh; ip++)
                {

                    newsintR_x[index] =
                        get_vel() * QMD_ddot (count, workR, incx, &beta_x[prjcount], incx);
                    newsintR_y[index] =
                        get_vel() * QMD_ddot (count, workR, incx, &beta_y[prjcount], incx);
                    newsintR_z[index] =
                        get_vel() * QMD_ddot (count, workR, incx, &beta_z[prjcount], incx);
#if !GAMMA_PT
                    newsintI_x[index] =
                        get_vel() * QMD_ddot (count, workI, incx, &beta_x[prjcount], incx);
                    newsintI_y[index] =
                        get_vel() * QMD_ddot (count, workI, incx, &beta_y[prjcount], incx);
                    newsintI_z[index] =
                        get_vel() * QMD_ddot (count, workI, incx, &beta_z[prjcount], incx);
#endif
                    prjcount += count;
                    index++;

                }               /*end for ip */
            }                   /*end if count */
            sta++;
        }                       /* end for istate */
    }

    if (count)
        delete[] workR;
    /*my_free(temp_psiR); */
# if !FDIFF_BETA
    if (size)
        delete[] beta_x;
#endif
}
