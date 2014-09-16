/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"
#include "common_prototypes.h"

void assign_weight (SPECIES * sp, int ion, fftw_complex * beptr, rmg_double_t * rtptr, rmg_double_t *Bweight)
{

    int idx, ix, iy, iz, *dvec;
    int idx1, docount;
    int *pidx, nldim;
    rmg_double_t *tem_array, *Btem_array;

    nldim = sp->nldim;
    idx = nldim * nldim * nldim;
    my_malloc (tem_array, idx, rmg_double_t);
    my_malloc (Btem_array, idx, rmg_double_t);
    weight_shift_center(sp, beptr);
    for(ix = 0; ix < nldim * nldim * nldim; ix++) 
        tem_array[ix] = creal(beptr[ix]);

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
                    rtptr[pidx[docount]] = creal(beptr[idx1]);
                    Bweight[pidx[docount]] = Btem_array[idx1];
                    if (cimag(beptr[idx1]) > 1.0e-8)
                    {
                        printf ("beptr[%d].im=%e\n", idx1, cimag(beptr[idx1]));
                        error_handler ("something wrong with the fourier transformation");
                    }
                    docount++;
                }

                idx++;
            }
        }
    }
    if (docount != pct.idxptrlen[ion])
    {
        printf ("docount = %d != %d = pct.idxptrlen[ion = %d]\n", docount, pct.idxptrlen[ion], ion);
        error_handler ("wrong numbers of projectors");
    }

    my_free(Btem_array);
    my_free(tem_array);
}
