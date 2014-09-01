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

//template void AssignWeight<double> (Kpoint<double> **, SPECIES *, int, fftw_complex *, double *, double *);
//template void AssignWeight<std::complex<double> >(Kpoint<std::complex<double>> **, SPECIES *, int, fftw_complex *, double *, double *);

void AssignWeight (Kpoint<double> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, double * rtptr, double *Bweight, double *Nlweight)
{

    int pbasis = kptr->pbasis;
    int nldim = sp->nldim;

    for(int idx = 0; idx < pbasis; idx++) rtptr[idx] = 0.0;
    if(pct.idxptrlen[ion] == 0) return;

    std::complex<double> *nbeptr = (std::complex<double> *)beptr;


    double *tem_array = new double[nldim * nldim * nldim];
    double *Btem_array = new double[nldim * nldim * nldim];

    for(int ix = 0; ix < nldim * nldim * nldim; ix++) 
        tem_array[ix] = std::real(nbeptr[ix]);

    app_cir_beta_driver (tem_array, Btem_array, nldim, nldim, 
            nldim, ct.kohn_sham_fd_order);

    int *pidx = pct.nlindex[ion];
    int *dvec = pct.idxflag[ion];
    int idx = 0;
    int docount = 0;
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

void AssignWeight (Kpoint<std::complex<double>> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, double * rtptr, std::complex<double> *Bweight, std::complex<double> *Nlweight)
{

    int pbasis = kptr->pbasis;
    int nldim = sp->nldim;

    for(int idx = 0; idx < pbasis; idx++) rtptr[idx] = 0.0;
    if(pct.idxptrlen[ion] == 0) return;

    std::complex<double> *nbeptr = (std::complex<double> *)beptr;


    double *tem_array = new double[nldim * nldim * nldim];
    double *Btem_array = new double[nldim * nldim * nldim];

    for(int ix = 0; ix < nldim * nldim * nldim; ix++) 
        tem_array[ix] = std::real(nbeptr[ix]);

    app_cir_beta_driver (tem_array, Btem_array, nldim, nldim, 
          nldim, ct.kohn_sham_fd_order);


    double *pR = pct.phaseptr[ion];
    pR += 2 * kptr->kidx * pbasis;
    double *pI = pR + pbasis;

    int *pidx = pct.nlindex[ion];
    int *dvec = pct.idxflag[ion];
    int idx = 0;
    int docount = 0;
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
                    Nlweight[pidx[docount]] = std::complex<double>(std::real(nbeptr[idx1]), 0.0);
                    Bweight[pidx[docount]] = std::complex<double>(Btem_array[idx1]*pR[pidx[docount]], Btem_array[idx1]*pI[pidx[docount]]);
//Bweight[pidx[docount]] = std::complex<double>(Btem_array[idx1], 0.0);
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
