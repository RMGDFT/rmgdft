/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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


template void GetDerweight<double> ( int ion, double *beta_x, double *beta_y, double *beta_z, 
        ION *iptr, fftw_plan p2, Kpoint<double> *kptr);
template void GetDerweight<std::complex<double> > ( int ion, std::complex<double> *beta_x, 
        std::complex<double> *beta_y, std::complex<double> *beta_z, ION *iptr, fftw_plan p2,
        Kpoint<std::complex<double>> *kptr);



template <typename OrbitalType> void GetDerweight (int ion, OrbitalType * beta_x, 
        OrbitalType * beta_y, OrbitalType * beta_z, ION * iptr, fftw_plan p2, 
        Kpoint<OrbitalType> *kptr)
{

    int ip, coarse_size, idx, nh;
    SPECIES *sp;
    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr_x, *fptr_y, *fptr_z;
    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> *beptr, *gbptr;
    int P0_BASIS;

    P0_BASIS = get_P0_BASIS();



    /* Get species type */
    sp = &ct.sp[iptr->species];

    /*Number of grid points on which fourier transform is done (in the corse grid) */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;

    nh = pct.idxptrlen[ion];


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    beptr = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * 2* coarse_size);
    if (beptr == NULL)
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    gbptr = beptr + coarse_size;



    /*Temporary pointer to the already calculated forward transform */
    fptr_x = (std::complex<double> *)sp->forward_derbeta_x;
    fptr_y = (std::complex<double> *)sp->forward_derbeta_y;
    fptr_z = (std::complex<double> *)sp->forward_derbeta_z;




    /*Calculate the phase factor */
    //find_phase (sp->nldim, iptr->nlcrds, iptr->fftw_phase_sin, iptr->fftw_phase_cos);

    /* Loop over radial projectors */
    for (ip = 0; ip < sp->num_projectors; ip++)
    {


        /*************** X ********************/
        /*Apply the phase factor */
        for (idx = 0; idx < coarse_size; idx++)
        {
            gbptr[idx] =
                (std::real(fptr_x[idx]) * iptr->fftw_phase_cos[idx] +
                 std::imag(fptr_x[idx]) * iptr->fftw_phase_sin[idx]) +
                (std::imag(fptr_x[idx]) * iptr->fftw_phase_cos[idx] -
                 std::real(fptr_x[idx]) * iptr->fftw_phase_sin[idx]) * I_t;
        }


        int idx1, idx2;

        /*Do the backwards transform */
        fftw_execute_dft (p2, reinterpret_cast<fftw_complex*>(gbptr), reinterpret_cast<fftw_complex*>(beptr));
        /*This takes and stores the part of beta that is useful for this PE */

        AssignDerweight (sp, ion, reinterpret_cast<fftw_complex*>(beptr), &beta_x[ip * P0_BASIS], kptr);


        /*************** Y ********************/
        /*Apply the phase factor */
        for (idx = 0; idx < coarse_size; idx++)
        {
            gbptr[idx] =
                (std::real(fptr_y[idx]) * iptr->fftw_phase_cos[idx] +
                 std::imag(fptr_y[idx]) * iptr->fftw_phase_sin[idx]) +
                (std::imag(fptr_y[idx]) * iptr->fftw_phase_cos[idx] -
                 std::real(fptr_y[idx]) * iptr->fftw_phase_sin[idx]) * I_t;
        }

        /*Do the backwards transform */
        fftw_execute_dft (p2, reinterpret_cast<fftw_complex*>(gbptr), reinterpret_cast<fftw_complex*>(beptr));
        /*This takes and stores the part of beta that is useful for this PE */
        AssignDerweight (sp, ion, reinterpret_cast<fftw_complex*>(beptr), &beta_y[ip * P0_BASIS], kptr);


        /*************** Z ********************/
        /*Apply the phase factor */
        for (idx = 0; idx < coarse_size; idx++)
        {
            gbptr[idx] =
                (std::real(fptr_z[idx]) * iptr->fftw_phase_cos[idx] +
                 std::imag(fptr_z[idx]) * iptr->fftw_phase_sin[idx]) +
                (std::imag(fptr_z[idx]) * iptr->fftw_phase_cos[idx] -
                 std::real(fptr_z[idx]) * iptr->fftw_phase_sin[idx]) * I_t;
        }

        /*Do the backwards transform */
        fftw_execute_dft (p2, reinterpret_cast<fftw_complex*>(gbptr), reinterpret_cast<fftw_complex*>(beptr));
        /*This takes and stores the part of beta that is useful for this PE */
        AssignDerweight (sp, ion, reinterpret_cast<fftw_complex*>(beptr), &beta_z[ip * P0_BASIS], kptr);



        /*Advance the temp pointers */
        fptr_x += coarse_size;
        fptr_y += coarse_size;
        fptr_z += coarse_size;

    }                           /*end for(ip = 0;ip < sp->num_projectors;ip++) */




    delete[] beptr;


}                               /* end get_weight */
