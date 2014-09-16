#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "FiniteDiff.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

template void AssignWeight<double> (Kpoint<double> *, SPECIES *, int, fftw_complex *, double *, double *, double *);
template void AssignWeight<std::complex<double> >(Kpoint<std::complex<double>> *, SPECIES *, int, fftw_complex *, double *, std::complex<double> *, std::complex<double> *);


template <typename KpointType>
void AssignWeight (Kpoint<KpointType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, double * rtptr, KpointType *Bweight, KpointType *Nlweight)
{

    Lattice *L = kptr->L;
    TradeImages *T = kptr->T;


    int pbasis = kptr->pbasis;
    int nldim = sp->nldim;
    KpointType ZERO_t(0.0);

    std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight;
    std::complex<double> *Bweight_C = (std::complex<double> *)Bweight;

    double *Nlweight_R = (double *)Nlweight;
    double *Bweight_R = (double *)Bweight;

    weight_shift_center(sp, beptr);

    for(int idx = 0; idx < pbasis; idx++) rtptr[idx] = 0.0;
    for(int idx = 0; idx < pbasis; idx++) Bweight[idx] = ZERO_t;
    for(int idx = 0; idx < pbasis; idx++) Nlweight[idx] = ZERO_t;
    if(pct.idxptrlen[ion] == 0) return;

    std::complex<double> *nbeptr = (std::complex<double> *)beptr;


    KpointType *tem_array = new KpointType[nldim * nldim * nldim]();
    KpointType *Btem_array = new KpointType[nldim * nldim * nldim]();
    std::complex<double> *tem_array_C = (std::complex<double> *)tem_array;
    std::complex<double> *Btem_array_C = (std::complex<double> *)Btem_array;
    double *Btem_array_R = (double *)Btem_array;


    for(int ix = 0; ix < nldim * nldim * nldim; ix++) {
        tem_array[ix] = std::real(nbeptr[ix]);
    }

    double *pR = pct.phaseptr[ion];
    pR += 2 * kptr->kidx * pbasis;
    double *pI = pR + pbasis;

    int *pidx = pct.nlindex[ion];
    int *dvec = pct.idxflag[ion];
    int idx = 0;
    int docount = 0;


    // Apply phase factor
    for (int ix = 0; ix < sp->nldim; ix++)
    {

        for (int iy = 0; iy < sp->nldim; iy++)
        {

            for (int iz = 0; iz < sp->nldim; iz++)
            {

                if (dvec[idx])
                {
                    int idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                    rtptr[pidx[docount]] = std::real(nbeptr[idx1]);
                    if(!ct.is_gamma) {
                        tem_array_C[idx1] = std::real(nbeptr[idx1]) * std::complex<double> (pR[pidx[docount]], -pI[pidx[docount]]);
                    }
                    docount++;
                }

                idx++;
            }
        }
    }
    if (docount != pct.idxptrlen[ion])
    {
        rmg_printf ("docount = %d != %d = pct.idxptrlen[ion = %d]\n", docount, pct.idxptrlen[ion], ion);
        rmg_error_handler (__FILE__, __LINE__, "wrong numbers of projectors");
    }


    // Apply B operator then map weights back
    AppCirDriverBeta (L, T, tem_array, Btem_array, nldim, nldim, nldim, ct.kohn_sham_fd_order);

    idx = 0;
    docount = 0;
    for (int ix = 0; ix < sp->nldim; ix++)
    {

        for (int iy = 0; iy < sp->nldim; iy++)
        {

            for (int iz = 0; iz < sp->nldim; iz++)
            {

                if (dvec[idx])
                {
                    int idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                    rtptr[pidx[docount]] = std::real(nbeptr[idx1]);
                    if(ct.is_gamma) {
                        Nlweight_R[pidx[docount]] = std::real(nbeptr[idx1]);
                        Bweight_R[pidx[docount]] = Btem_array_R[idx1];
                    }
                    else {
                        Nlweight_C[pidx[docount]] = std::complex<double>(std::real(nbeptr[idx1]) * pR[pidx[docount]], -std::real(nbeptr[idx1]) * pI[pidx[docount]]);
                        Bweight_C[pidx[docount]] = Btem_array_C[idx1];
                    }
                    docount++;
                }

                idx++;
            }
        }
    }
    if (docount != pct.idxptrlen[ion])
    {
        rmg_printf ("docount = %d != %d = pct.idxptrlen[ion = %d]\n", docount, pct.idxptrlen[ion], ion);
        rmg_error_handler (__FILE__, __LINE__, "wrong numbers of projectors");
    }


    delete [] Btem_array;
    delete [] tem_array;



    
}
