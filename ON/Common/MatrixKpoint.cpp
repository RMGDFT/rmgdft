
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex>

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_control.h"
#include "transition.h"
#include "blacs.h"

void MatrixKpointPhase (STATE *states, int *desca, std::vector<int> &min_index)
{

    int mycol, myrow, nprow, npcol;
    int ictxt=desca[1], m = desca[2], n = desca[3], mb=desca[4], nb=desca[5], mxllda = desca[8];

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int izero = 0;
    int mxli = numroc(&m, &mb, &myrow, &izero, &nprow);
    int mxlocc = numroc(&n, &nb, &mycol, &izero, &npcol);


    double crds[3], xtal[3], xtal_t[3];
    std::vector<double> dists;
    dists.resize(27);
    std::complex<double> I(0.0, 1.0);


    for(int j =0; j < mxli; j++)
    {
        for(int k=0; k < mxlocc; k++)
        {

            int jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb;
            int kk = (k/nb ) * npcol * nb + mycol * nb + k - k/nb * nb;

            crds[0] = states[jj].crds[0] - states[kk].crds[0];
            crds[1] = states[jj].crds[1] - states[kk].crds[1];
            crds[2] = states[jj].crds[2] - states[kk].crds[2];

            Rmg_L.to_crystal_vector(xtal, crds);

            int idx = 0;
            for(int ix = -1; ix <= 1; ix++)
            {
                for(int iy = -1; iy <= 1; iy++)
                {
                    for(int iz = -1; iz <= 1; iz++)
                    {
                        xtal_t[0] = xtal[0] + ix;
                        xtal_t[1] = xtal[1] + iy;
                        xtal_t[2] = xtal[2] + iz;

                        dists[idx] = Rmg_L.metric(xtal_t) ;
                        idx++;
                    }
                }
            }

            std::vector<double>::iterator min_iter = std::min_element(dists.begin(), dists.end());
            min_index[k *mxllda + j] = std::distance(dists.begin(), min_iter);
        }
    }
}

void MatrixKpoint (STATE *states, double *Hij, double *Sij, int *desca, 
        std::complex<double> *Hk, std::complex<double> *Sk, double kpt[3])
{

    int mycol, myrow, nprow, npcol;
    int ictxt=desca[1], m = desca[2], n = desca[3], mb=desca[4], nb=desca[5], mxllda = desca[8];

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int izero = 0;
    int mxli = numroc(&m, &mb, &myrow, &izero, &nprow);
    int mxlocc = numroc(&n, &nb, &mycol, &izero, &npcol);

    static std::vector<int> min_index;
    static std::complex<double> phase_k[27];


    double crds[3], xtal[3], xtal_t[3];
    std::vector<double> dists;
    dists.resize(27);
    std::complex<double> I(0.0, 1.0);


    if(min_index.size() == 0)
    {
        min_index.resize(mxllda * mxlocc);
        int idx = 0;
        for(int ix = -1; ix <= 1; ix++)
        {
            for(int iy = -1; iy <= 1; iy++)
            {
                for(int iz = -1; iz <= 1; iz++)
                {
                    phase_k[idx] = std::exp(+I * (ix * kpt[0] + iy * kpt[1] + iz * kpt[2]));
                    idx++;
                }
            }
        }


        for(int j =0; j < mxli; j++)
        {
            for(int k=0; k < mxlocc; k++)
            {

                int jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb;
                int kk = (k/nb ) * npcol * nb + mycol * nb + k - k/nb * nb;

                crds[0] = states[jj].crds[0] - states[kk].crds[0];
                crds[1] = states[jj].crds[1] - states[kk].crds[1];
                crds[2] = states[jj].crds[2] - states[kk].crds[2];

                Rmg_L.to_crystal_vector(xtal, crds);

                int idx = 0;
                for(int ix = -1; ix <= 1; ix++)
                {
                    for(int iy = -1; iy <= 1; iy++)
                    {
                        for(int iz = -1; iz <= 1; iz++)
                        {
                            xtal_t[0] = xtal[0] + ix;
                            xtal_t[1] = xtal[1] + iy;
                            xtal_t[2] = xtal[2] + iz;

                            dists[idx] = Rmg_L.metric(xtal_t) ;
                            idx++;
                        }
                    }
                }

                std::vector<double>::iterator min_iter = std::min_element(dists.begin(), dists.end());
                min_index[k *mxllda + j] = std::distance(dists.begin(), min_iter);
            }
        }
    }

    for(int idx = 0; idx < mxllda * mxlocc; idx++)
    {
        Hk[idx] = Hij[idx] * phase_k[ min_index[idx] ];
        Sk[idx] = Sij[idx] * phase_k[ min_index[idx] ];
    }

}

