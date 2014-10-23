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
    OrbitalType *beta_x, *beta_y, *beta_z, temx, temy, temz;
    double *temp_psiR;
    SPECIES *sp;
    int P0_BASIS;
    double *workI, *pI, *temp_psiI, *pR;

    nh = ct.sp[ct.ions[ion].species].nh;
    sp = &ct.sp[ct.ions[ion].species];
    P0_BASIS = get_P0_BASIS();
    int xoff = get_PX_OFFSET();

    count = pct.idxptrlen[ion];
    pidx = pct.nlindex[ion];

    double tem, tem1;
    if (!count) return;

    workR = new double[2* P0_BASIS];
    workI = workR + P0_BASIS;


    size = nh * P0_BASIS;

    beta_x = new OrbitalType[ 3 * size]; 

    beta_y = beta_x + size;
    beta_z = beta_y + size;



    for (kidx = 0; kidx < ct.num_kpts; kidx++)
    {

        for (idx = 0; idx < size; idx++)
        {
            beta_x[idx] = 0.0;
            beta_y[idx] = 0.0;
            beta_z[idx] = 0.0;
        }

        /*partial_beta(ion, beta_x, beta_y,beta_z, iptr, p1, p2); */
        GetDerweight (ion, beta_x, beta_y, beta_z, iptr, p2, Kptr[kidx]);

        for (istate = 0; istate < ct.num_states; istate++)
        {


            /*Gather_psi is not necessary, getting pointer should be enough
              gather_psi(temp_psiR, temp_psiI, sta, 0); */

            index =
                ion * ct.num_kpts * ct.num_states * ct.max_nl +
                kidx * ct.num_states * ct.max_nl + istate * ct.max_nl;
            for (ip = 0; ip < nh; ip++)
            {

                temx = 0.0; temy = 0.0; temz = 0.0;
            for (idx = 0; idx < P0_BASIS; idx++)
            {
                temx += Kptr[kidx]->Kstates[istate].psi[idx] * std::conj(beta_x[ip * P0_BASIS + idx]);
                temy += Kptr[kidx]->Kstates[istate].psi[idx] * std::conj(beta_y[ip * P0_BASIS + idx]);
                temz += Kptr[kidx]->Kstates[istate].psi[idx] * std::conj(beta_z[ip * P0_BASIS + idx]);
            }

                newsintR_x[index+ip] = get_vel() * std::real(temx);
                newsintR_y[index+ip] = get_vel() * std::real(temy);
                newsintR_z[index+ip] = get_vel() * std::real(temz);
                if(!ct.is_gamma)
                {
                    newsintI_x[index+ip] =get_vel() * std::imag(temx);
                    newsintI_y[index+ip] =get_vel() * std::imag(temy);
                    newsintI_z[index+ip] =get_vel() * std::imag(temz);
                }

            }               /*end for ip */
        }                   /*end if count */
    }                       /* end for istate */

    delete[] workR;
    delete[] beta_x;
}
