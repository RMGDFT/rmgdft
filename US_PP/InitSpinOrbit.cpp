#include <stdio.h>
#include <float.h>
#include <cstdlib>
#include <math.h>
#include "blas.h"
#include "transition.h"
#include "const.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>

static void InitUmm(int lmax, std::complex<double> *Umm);
static void Init_fcoef(SPECIES &sp, std::complex<double> *Umm, int tot_LM);
void InitSpinOrbit ()
{
    int tot_LM = (ct.max_l +1) *(ct.max_l +1);
    std::complex<double> *Umm = new std::complex<double>[tot_LM * tot_LM]();
    InitUmm(ct.max_l, Umm);


    for(int isp = 0; isp < ct.num_species; isp++)
    {
        SPECIES *sp = &Species[isp];
        Init_fcoef(*sp, Umm, tot_LM);
    }
    delete [] Umm;
}

static void Init_fcoef(SPECIES &sp, std::complex<double> *Umm, int tot_LM)
{
    for(int ih = 0; ih < sp.nh; ih++)
    for(int jh = 0; jh < sp.nh; jh++)
    {
        sp.fcoef_so[ih][jh][0] = 0.0;
        sp.fcoef_so[ih][jh][1] = 0.0;
        sp.fcoef_so[ih][jh][2] = 0.0;
        sp.fcoef_so[ih][jh][3] = 0.0;
    }
    
    if(!sp.is_spinorb)  return;

    double alpha_up, alpha_dn;
    std::complex<double> Umi_up, Umj_up, Umi_dn, Umj_dn;
    for(int ih = 0; ih < sp.nh; ih++)
    {
        int li = sp.nhtol[ih];
        double ji = sp.nhtoj[ih];
        int lmi = sp.nh_l2m[ih];

        for(int jh = 0; jh < sp.nh; jh++)
        {
            int lj = sp.nhtol[jh];
            double jj = sp.nhtoj[jh];
            int lmj = sp.nh_l2m[jh];

            //  fcoef != 0.0 only when l and j are the same. 
            if( li != lj || std::abs(ji - jj) > 1.0-5) continue;
         // case: j = l + 1/2, mj =[-j,j], m = mj-1/2 = [-l-1, l]
            if( std::abs(ji -li - 0.5) < 1.0e-5)
            {
                for(int m = -li -1; m <= li; m++)
                {
                    alpha_up = std::sqrt( (li + m + 1.0)/(2*li + 1.0));
                    alpha_dn = std::sqrt( (li - m )/(2*li + 1.0));
                    int lmm = li * li + li-m;
                    if(m  < -li) 
                    {
                        Umi_up = 0.0;
                        Umj_up = 0.0;
                    }
                    else
                    {
                        Umi_up = Umm[lmm * tot_LM + lmi];
                        Umj_up = std::conj(Umm[lmm * tot_LM + lmj]);
                    }

                    lmm = li * li + li - (m+1);
                    if( m +1 > li) 
                    {
                        Umi_dn = 0.0;
                        Umj_dn = 0.0;
                    }
                    else
                    {
                        Umi_dn = Umm[lmm * tot_LM + lmi];
                        Umj_dn = std::conj(Umm[lmm * tot_LM + lmj]);
                    }

                    sp.fcoef_so[ih][jh][0] += alpha_up *alpha_up * Umi_up * Umj_up; 
                    sp.fcoef_so[ih][jh][1] += alpha_up *alpha_dn * Umi_up * Umj_dn; 
                    sp.fcoef_so[ih][jh][2] += alpha_dn *alpha_up * Umi_dn * Umj_up; 
                    sp.fcoef_so[ih][jh][3] += alpha_dn *alpha_dn * Umi_dn * Umj_dn; 

                }
            }
        //  case: j = l - 1/2, mj =[-j,j], m = mj+1/2 = [-l+1, l]
            else if( std::abs(ji -li + 0.5) < 1.0e-5)
            {
                for(int m = -li+1; m <= li; m++)
                {
                    alpha_up = std::sqrt( (li - m + 1.0)/(2*li + 1.0));
                    alpha_dn = -std::sqrt( (li + m )/(2*li + 1.0));
                    int lmm = li * li + li-(m-1);
                    Umi_up = Umm[lmm * tot_LM + lmi];
                    Umj_up = std::conj(Umm[lmm * tot_LM + lmj]);
                    lmm = li * li + li-m;
                    Umi_dn = Umm[lmm * tot_LM + lmi];
                    Umj_dn = std::conj(Umm[lmm * tot_LM + lmj]);

                    sp.fcoef_so[ih][jh][0] += alpha_up *alpha_up * Umi_up * Umj_up; 
                    sp.fcoef_so[ih][jh][1] += alpha_up *alpha_dn * Umi_up * Umj_dn; 
                    sp.fcoef_so[ih][jh][2] += alpha_dn *alpha_up * Umi_dn * Umj_up; 
                    sp.fcoef_so[ih][jh][3] += alpha_dn *alpha_dn * Umi_dn * Umj_dn; 

                }
            }

        }
    }

}
static void InitUmm(int lmax, std::complex<double> *Umm)
{
    // calculate transfor matrix from real  spherical harmonic to complex harmonics

    // real harmonic are from Ylm.cpp with index L * L + M, M = [0, 2L+1)
    // complex harmoics: index L*L + M (if M < 0, M = M + 2 *L+1)

    // total number of harmonic for LM 
    int tot_LM = (lmax +1) *(lmax +1);
    int idx, idx1, lm, num_lm;
    double r[3];

    double *ylm_real = new double[tot_LM * tot_LM];
    double *ylm_invert = new double[tot_LM * tot_LM];
    std::complex<double> *ylm_complex = new std::complex<double>[tot_LM * tot_LM];


    std::srand(224);
    //  generating tot_LM random r(3) and for each r, calculating Ylm(r) for all LMs
    for (int i = 0; i < tot_LM; i++)
    {

        double theta = std::rand()/(double)RAND_MAX * PI;
        double phi = std::rand()/(double)RAND_MAX * twoPI;


        r[0] = sin(theta) * cos(phi);
        r[1] = sin(theta) * sin(phi);
        r[2] = cos(theta);

        for(int L = 0; L <= lmax; L++)
            for(int M = 0; M < 2*L+1; M++)
            {
                lm = L * L + M;

                ylm_real[i * tot_LM + lm] = Ylm(L, M, r);
                ylm_complex[i * tot_LM + lm] = boost::math::spherical_harmonic(L, L-M, theta, phi);
            }
    }

    for(int i = 0; i < tot_LM * tot_LM; i++) ylm_invert[i] = 0.0;
    for(int i = 0; i < tot_LM; i++) ylm_invert[i * tot_LM + i] = 1.0;

    int info, ipvt[tot_LM];
    dgesv (&tot_LM, &tot_LM, ylm_real, &tot_LM, ipvt, ylm_invert, &tot_LM, &info);

    for(int i = 0; i < tot_LM; i++)
        for(int j = 0; j < tot_LM; j++)
            for(int k = 0; k < tot_LM; k++)
                Umm[i * tot_LM + j] += ylm_complex[k * tot_LM + j] * ylm_invert[i * tot_LM + k];

    if(pct.gridpe == 0 && 0)
    {
        for(int i = 0; i < tot_LM; i++)
        {
            printf("\n real U: ");
            for(int j = 0; j < tot_LM; j++)
                printf(" %e ", std::real(Umm[i*tot_LM + j]));
        }

        for(int i = 0; i < tot_LM; i++)
        {
            printf("\n imag U: ");
            for(int j = 0; j < tot_LM; j++)
                printf(" %e ", std::imag(Umm[i*tot_LM + j]));
        }
    }

    delete [] ylm_real;
    delete [] ylm_invert;
    delete [] ylm_complex;

}


