#include <stdio.h>
#include <float.h>
#include <cstdlib>
#include <math.h>
#include "blas.h"
#include "transition.h"

void InitClebschGordan (int lmax, double *ap, int *lpx, int *lpl)
{
   

    // calculate the expansion coefficients for real harmonic, similar to ClebschGordan for spherical harmonic
    // we follow the method in pwscf 
    //  Y_lmi(r) * Y_lmj(r) = Sum_LM ap[LM][lmi][lmj] *Y_LM(r)
  

    // total number of harmonic for LM 
    int tot_LM = (2 * lmax +1) *(2 * lmax +1);
    int idx, idx1, lm, num_lm;
    double r[3];

    double *ylm_array = new double[tot_LM * tot_LM];
    double *ylm_invert = new double[tot_LM * tot_LM];
    double *ylm_tem = new double[tot_LM * tot_LM];

    num_lm = (lmax+1) * (lmax+1);
    for(int i = 0; i < num_lm *num_lm * tot_LM; i++) 
        lpl[i] = 0;



    std::srand(224);
    //  generating tot_LM random r(3) and for each r, calculating Ylm(r) for all LMs
    for (int i = 0; i < tot_LM; i++)
    {

        r[0] = (2.0 *std::rand())/RAND_MAX -1.0;
        r[1] = (2.0 *std::rand())/RAND_MAX -1.0;
        r[2] = (2.0 *std::rand())/RAND_MAX -1.0;

        for(int L = 0; L <= lmax *2; L++)
            for(int M = 0; M < 2*L+1; M++)
            {
                lm = L * L + M;

                ylm_array[i * tot_LM + lm] = Ylm(L, M, r);
                ylm_tem[i * tot_LM + lm] = ylm_array[i * tot_LM + lm];
            }
    }

    for(int i = 0; i < tot_LM * tot_LM; i++) ylm_invert[i] = 0.0;
    for(int i = 0; i < tot_LM; i++) ylm_invert[i * tot_LM + i] = 1.0;

    int info, ipvt[tot_LM];
    dgesv (&tot_LM, &tot_LM, ylm_tem, &tot_LM, ipvt, ylm_invert, &tot_LM, &info);


    for(int il = 0; il < num_lm; il++)
        for(int ik = 0; ik < num_lm; ik++)
        {
            lpx[il*num_lm + ik] = 0;
            for(int ilp = 0; ilp < tot_LM; ilp++)
            {
                idx = ilp * num_lm *num_lm + il*num_lm + ik;
                ap[idx] = 0.0;
                for(int ir = 0; ir < tot_LM; ir++)
                    ap[idx] += ylm_invert[ilp * tot_LM + ir] *ylm_array[ir * tot_LM + il] * ylm_array[ir * tot_LM + ik];
                if(std::abs(ap[idx]) > 1.0e-3)
                {  
                    idx1 = il*num_lm + ik;
                    lpl[idx1 * tot_LM + lpx[idx1]] = ilp;
                    lpx[idx1] += 1;
                }
            }

        }


    delete [] ylm_array;
    delete [] ylm_invert;
    delete [] ylm_tem;
}
