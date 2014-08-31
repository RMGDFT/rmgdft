#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

template void AssignWeight<double> (Kpoint<double> **, SPECIES *, int, fftw_complex *, double *, double *);
template void AssignWeight<std::complex<double> >(Kpoint<std::complex<double>> **, SPECIES *, int, fftw_complex *, double *, double *);

template <typename KpointType>
void AssignWeight (Kpoint<KpointType> **Kptr, SPECIES * sp, int ion, fftw_complex * beptr, double * rtptr, double *Bweight)
{

    int idx, ix, iy, iz, *dvec;
    int idx1, docount;
    int *pidx, nldim;
    double *tem_array, *Btem_array;
    std::complex<double> *nbeptr = (std::complex<double> *)beptr;

    nldim = sp->nldim;
    idx = nldim * nldim * nldim;
    tem_array = new double[idx];
    Btem_array = new double[idx];
    for(ix = 0; ix < nldim * nldim * nldim; ix++) 
        tem_array[ix] = std::real(nbeptr[ix]);

    app_cir_beta_driver (tem_array, Btem_array, nldim, nldim, 
            nldim, ct.kohn_sham_fd_order);

    for(idx = 0; idx < get_P0_BASIS(); idx++) rtptr[idx] = 0.0;
    if(pct.idxptrlen[ion] == 0) return;
    pidx = pct.nlindex[ion];
    dvec = pct.idxflag[ion];
    idx = docount = 0;
    for (ix = 0; ix < sp->nldim; ix++)
    {

        for (iy = 0; iy < sp->nldim; iy++)
        {

            for (iz = 0; iz < sp->nldim; iz++)
            {

                if (dvec[idx])
                {
                    idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                    rtptr[pidx[docount]] = std::real(nbeptr[idx1]);
                    Bweight[pidx[docount]] = Btem_array[idx1];
                    if (std::imag(nbeptr[idx1]) > 1.0e-8)
                    {
                        rmg_printf ("beptr[%d].im=%e\n", idx1, std::imag(nbeptr[idx1]));
                        rmg_error_handler (__FILE__, __LINE__, "something wrong with the fourier transformation");
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
