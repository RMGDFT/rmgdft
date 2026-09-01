/************************** SVN Revision Information **************************
 **    $Id: get_nlop.c 3560 2016-05-10 14:27:30Z ebriggs $    **
 ******************************************************************************/

/*

   Sets up the ket part of the non-local operators.


 */


//#include <sys/types.h>
//#include <sys/stat.h>

#include <fcntl.h>
//#include <libgen.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex>
#include <fftw3.h>
#include <sys/mman.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>


#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"
#include "RmgTimer.h"
#include "RmgException.h"


void GetNlop_on(void)
{
    int ion, idx, ip;
    size_t prjcount;
    double *beta;
    SPECIES *sp;
    ION *iptr;
    fftw_plan p2;
    int coarse_size;
    double vect[3], nlcrds[3];

    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;
    std::complex<double> *beptr, *gbptr;

    /*Do forward transform for each species and store results on the coarse grid */
    if(ct.localize_projectors)
    {
        for(auto& sp : Species) sp.InitWeights (ct.localize_projectors);
    }
    else
    {
        throw RmgFatalException() << "delocalized projector not supported" << __FILE__ << " at line " << __LINE__ << "\n";
    }

    /*The same for derivative of beta */
    //init_derweight ();


    /*Get memory to store the phase factor applied to the forward Fourier transform
     *      * and to store the backwards transform*/
    beptr = new std::complex<double>[2 * ct.max_nlpoints];
    gbptr = beptr + ct.max_nlpoints;

    std::complex<double> *fftw_phase = new std::complex<double>[ct.max_nlpoints]; 


    /*
     * PROJECTOR_SPACE = ct.max_nlpoints * ct.max_nl;
     */

    MPI_Barrier(pct.img_comm);

    /* Loop over all the ions on this processor */

    mkdir("PROJECTORS",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    // beta = projectors;

    for (ion = pct.imgpe; ion < ct.num_ions; ion+=pct.image_npes[pct.thisimg])
    {

        /* Generate ion pointer */
        iptr = &Atoms[ion];

        /* Get species type */
        sp = &Species[iptr->species];
        beta = new double[sp->num_projectors * ct.max_nlpoints];

        p2 = fftw_plan_dft_3d(sp->nldim, sp->nldim, sp->nldim, reinterpret_cast<fftw_complex*>(gbptr), 
                reinterpret_cast<fftw_complex*>(beptr), FFTW_BACKWARD, FFTW_ESTIMATE);

        /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
        /*xtal vector between ion and left bottom corner of the box */

        vect[0] = iptr->xtal[0] - iptr->nlxcstart;
        vect[1] = iptr->xtal[1] - iptr->nlycstart;
        vect[2] = iptr->xtal[2] - iptr->nlzcstart;

        /*Substract vector between left bottom corner of the box and center of the box */
        vect[0] -= (sp->nldim / 2) / (double) get_NX_GRID();
        vect[1] -= (sp->nldim / 2) / (double) get_NY_GRID();
        vect[2] -= (sp->nldim / 2) / (double) get_NZ_GRID();

        /*The vector we are looking for should be */
        to_cartesian (vect, nlcrds);
        coarse_size = sp->nldim *sp->nldim *sp->nldim ;

        /*Calculate the phase factor */
        FindPhase (sp, sp->nldim, sp->nldim, sp->nldim, nlcrds, fftw_phase);

        /*Temporary pointer to the already calculated forward transform */
        fptr = (std::complex<double> *)sp->forward_beta;

        /* Loop over radial projectors */
        prjcount = 0;
        for (ip = 0; ip < sp->num_projectors; ip++)
        {

            /*Apply the phase factor   */
            for (idx = 0; idx < coarse_size; idx++)
            {
                gbptr[idx] = fptr[idx] * std::conj(fftw_phase[idx]);
            }


            /*Do the backwards transform */
            fftw_execute_dft (p2, reinterpret_cast<fftw_complex*>(gbptr), reinterpret_cast<fftw_complex*>(beptr));
            /*This takes and stores the part of beta that is useful for this PE */
            assign_weight_on (sp, reinterpret_cast<fftw_complex*>(beptr), &beta[prjcount]);



            fptr += coarse_size;
            prjcount += ct.max_nlpoints;

        }                       /* end for ip */

        fftw_destroy_plan(p2);
        std::string newname;
        newname = std::string("PROJECTORS/ion_") + std::to_string(ion);
        int amode = S_IREAD | S_IWRITE;
        int fhand = open(newname.c_str(), O_CREAT | O_TRUNC | O_RDWR, amode);
        if (fhand < 0) 
            rmg::error("error open file");

        ssize_t size = (ssize_t)sp->num_projectors * (ssize_t)ct.max_nlpoints * sizeof(double);
        rmg::writefile(fhand, beta, size);
        close(fhand);

        delete [] beta;

    }                           /* end for ion */
    MPI_Barrier(pct.img_comm);

    delete [] beptr;
    delete [] fftw_phase;
    // Must fix this EMIL
    //


#if	DEBUG
    rmg::printlog("PE: %d leave  get_nlop ...\n", pct.gridpe);
    fflush(NULL);
#endif

    if (pct.gridpe == 0)
    {

        rmg::printlog(" get_nlop.c  done\n");

    }                           /* end if */
    /*    MPI_Barrier(pct.img_comm); */
    fflush(NULL);

}                               /* end get_nlop */
