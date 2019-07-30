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

static void init_alloc_nonloc_mem (void);

void GetNlop_on(void)
{
    int ion, idx, ip;
    int tot_prj, index;
    size_t PROJECTOR_SPACE;
    size_t prjcount;
    double *beta;
    SPECIES *sp;
    ION *iptr;
    fftw_plan p2;
    int overlap;
    int coarse_size, st1;
    double vect[3], nlcrds[3];

    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;
    std::complex<double> *beptr, *gbptr;

    init_alloc_nonloc_mem ();


    /*Do forward transform for each species and store results on the coarse grid */
    InitLocalizedWeight ();
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


           

    /*  get total number of projectors on this processor */
    /*  pct.n_ion_center: number of ions whose nl projector overlap
     *  with the states on this processor */

    pct.n_ion_center = 0;
    tot_prj = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        overlap = 0;
        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        {
            index = (st1 - ct.state_begin) * ct.num_ions + ion;
            if (ion_orbit_overlap_region_nl[index].flag == 1)
                overlap = 1;
        }
        if (overlap == 1)
        {
            pct.ionidx[pct.n_ion_center] = ion;
            pct.n_ion_center += 1;
            tot_prj += ct.sp[Atoms[ion].species].num_projectors;
        }
    }

    PROJECTOR_SPACE = (size_t)ct.max_nlpoints * (size_t)tot_prj;


//    printf("\n proj  %d %d %lu\n", ct.max_nlpoints, tot_prj, PROJECTOR_SPACE);
    std::string newpath;

    if(ct.nvme_weights)
    {
        if(ct.nvme_weight_fd != -1) close(ct.nvme_weight_fd);

        newpath = ct.nvme_weights_path + std::string("rmg_weight") + std::to_string(pct.spinpe) +
                  std::to_string(pct.kstart) + std::to_string(pct.gridpe);
        ct.nvme_weight_fd = FileOpenAndCreate(newpath, O_RDWR|O_CREAT|O_TRUNC, (mode_t)0600);

        projectors = (double *)CreateMmapArray(ct.nvme_weight_fd, PROJECTOR_SPACE*sizeof(double));
        if(!projectors) rmg_error_handler(__FILE__,__LINE__,"Error: CreateMmapArray failed for GetNlop_on. \n");
        madvise(projectors, PROJECTOR_SPACE*sizeof(double), MADV_SEQUENTIAL);
        
    }
    else
    {
        if (projectors != NULL)
            delete []projectors;
        projectors = new double[PROJECTOR_SPACE];
    }

    /*allocate memorry for weight factor of partial_beta/partial_x */


    for (ion = 0; ion < ct.num_ions; ion++)
        pct.prj_per_ion[ion] = ct.sp[Atoms[ion].species].num_projectors;

    /* Loop over all the ions on this processor */

    mkdir("PROJECTORS",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    beta = projectors;

    for (ion = pct.gridpe; ion < ct.num_ions; ion+=pct.grid_npes)
    {

        /* Generate ion pointer */
        iptr = &Atoms[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        //        fftw_import_wisdom_from_string(sp->backward_wisdom);
        p2 = fftw_plan_dft_3d(sp->nldim, sp->nldim, sp->nldim, reinterpret_cast<fftw_complex*>(gbptr), 
                reinterpret_cast<fftw_complex*>(beptr), FFTW_BACKWARD, FFTW_ESTIMATE);
        //        fftw_forget_wisdom();

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
        FindPhase (sp->nldim, sp->nldim, sp->nldim, nlcrds, fftw_phase);

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
            rmg_error_handler (__FILE__, __LINE__,"error open file");

        ssize_t size = (ssize_t)sp->num_projectors * (ssize_t)ct.max_nlpoints * sizeof(double);
        write(fhand, beta, size);
        close(fhand);


    }                           /* end for ion */

    delete [] beptr;
    delete [] fftw_phase;
    // Must fix this EMIL
    //
    
    MPI_Barrier(pct.grid_comm);
    prjcount = 0;
    for (unsigned int ion1 = 0; ion1 < pct.n_ion_center; ion1++)
    {
        ion = pct.ionidx[ion1];
        /* Generate ion pointer */
        iptr = &Atoms[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];
        std::string newname;
        newname = std::string("PROJECTORS/ion_") + std::to_string(ion);
        int fhand = open(newname.c_str(), O_RDWR, S_IREAD | S_IWRITE);
        ssize_t size = (ssize_t)sp->num_projectors * (ssize_t)ct.max_nlpoints * sizeof(double);
        ssize_t read_size = read(fhand, &beta[prjcount], size);
        if(read_size != size)
            rmg_error_handler (__FILE__, __LINE__,"error reading");
        prjcount += sp->num_projectors * ct.max_nlpoints;
    }

    MPI_Barrier(pct.grid_comm);
    

#if	DEBUG
    printf("PE: %d leave  get_nlop ...\n", pct.gridpe);
    fflush(NULL);
#endif

    if (pct.gridpe == 0)
    {

        printf(" get_nlop.c  done\n");

    }                           /* end if */
    /*    MPI_Barrier(pct.img_comm); */
    fflush(NULL);

}                               /* end get_nlop */


static void init_alloc_nonloc_mem (void)
{
    int ion;


    pct.Qindex = new int *[ct.num_ions];
    pct.Qdvec = new int *[ct.num_ions];
    pct.Qidxptrlen = new int[ct.num_ions];
    pct.lptrlen = new int[ct.num_ions];
    pct.ionidx = new int[ct.num_ions];
    pct.prj_per_ion = new int[ct.num_ions];

    pct.augfunc = new double *[ct.num_ions];
    pct.dnmI = new double *[ct.num_ions];
    pct.qqq = new double *[ct.num_ions];


    pct.ionidx_loc = new int[ct.num_ions];
    pct.prj_ptr = new int[ct.num_ions];

    for (ion = 0; ion < ct.num_ions; ion++)
    {


        pct.Qidxptrlen[ion] = 0;
        pct.lptrlen[ion] = 0;
        pct.ionidx[ion] = 0;
        pct.ionidx_loc[ion] = 0;
        pct.prj_per_ion[ion] = 0;

        pct.Qindex[ion] = NULL;
        pct.Qdvec[ion] = NULL;

        pct.augfunc[ion] = NULL;
        pct.dnmI[ion] = NULL;
        pct.qqq[ion] = NULL;

    }                           /*end for(ion=0; ion<ct.num_ions; ion++) */

}                               /*end init_alloc_nonloc_mem */


