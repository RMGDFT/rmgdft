/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

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
#include "RmgException.h"
#include "RmgShm.h"


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void InitWeight (void)
{

    int ip, prjcount, isp, size, tot_proj;
    SPECIES *sp;
    fftw_plan p1;
    fftw_complex *in, *out;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;

    // get tot number of projectors and their information

    tot_proj = 0;
    for (isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        sp = &ct.sp[isp];


        /*Loop over all betas to calculate num of projectors for given species */
        prjcount = 0;
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            for(int m = 0; m < 2*sp->llbeta[ip]+1; m++)
            {
                proj.species = isp;
                proj.ip = ip;
                proj.l = sp->llbeta[ip];
                proj.m = m;
                proj.proj_index = prjcount;
                proj_iter.push_back(proj);
                prjcount++;
            }
        }

        size = sp->nldim * sp->nldim * sp->nldim;
        if(!ct.localize_projectors) size = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
        sp->num_projectors = prjcount;

        /*This array will store forward fourier transform on the coarse grid for all betas */
        if(pct.procs_per_host > 1) {
            char sname[256];
            snprintf(sname, sizeof(sname), "RMG_ForwardBeta_%s", sp->atomic_symbol);
            sp->forward_beta = (fftw_complex *)AllocSharedMemorySegment(sname, sizeof(fftw_complex) * sp->num_projectors * size);
        }
        if(!sp->forward_beta) {
            sp->forward_beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_projectors * size);
        }

        if (sp->forward_beta == NULL)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

        tot_proj += prjcount;
    }


    for(int iproj = pct.gridpe; iproj < tot_proj; iproj+=pct.grid_npes)
    {
        proj = proj_iter[iproj];

        sp = &ct.sp[proj.species];
 
        size = sp->nldim * sp->nldim * sp->nldim;
        if(!ct.localize_projectors) size = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

        
        /*This is something we need to do only once per species, so do not use wisdom */
        if(ct.localize_projectors) {
            in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
            out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
        }
        else {
            in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ct.nxfgrid * ct.nxfgrid * ct.nxfgrid * get_NX_GRID() *  get_NY_GRID() * get_NZ_GRID());
            out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ct.nxfgrid * ct.nxfgrid * ct.nxfgrid * get_NX_GRID() *  get_NY_GRID() * get_NZ_GRID());
        }


        if(!in || !out)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

        if(ct.localize_projectors) {
            p1 = fftw_plan_dft_3d (sp->nlfdim, sp->nlfdim, sp->nlfdim, in, out, FFTW_FORWARD, FFTW_MEASURE);
        }
        else {
            p1 = fftw_plan_dft_3d (get_NX_GRID()*ct.nxfgrid, get_NY_GRID()*ct.nxfgrid, get_NZ_GRID()*ct.nxfgrid, in, out, FFTW_FORWARD, FFTW_MEASURE);
        }

        InitWeightOne(sp, &sp->forward_beta[proj.proj_index * size], proj.ip, proj.l, proj.m, p1);

        fftw_destroy_plan (p1);
        fftw_free(out);
        fftw_free(in);

    }                           /* end for */

    int root;
    for(int iproj = 0; iproj < tot_proj; iproj++)
    {
        proj = proj_iter[iproj];
        sp = &ct.sp[proj.species];
        root = iproj % pct.grid_npes;
        size = sp->nldim * sp->nldim * sp->nldim;
        if(!ct.localize_projectors) size = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

        MPI_Bcast(&sp->forward_beta[proj.proj_index * size], 2*size, MPI_DOUBLE, root, pct.grid_comm);
    }
        

    // Next if using non-localized projectors we need to remap the global forward beta into domain-decomposed
    // forward beta
    if(ct.localize_projectors) return;


    int ilo = get_PX_OFFSET();
    int jlo = get_PY_OFFSET();
    int klo = get_PZ_OFFSET();
    int ihi = ilo + get_PX0_GRID();
    int jhi = jlo + get_PY0_GRID();
    int khi = klo + get_PZ0_GRID();
    int dimx = get_NX_GRID();
    int dimy = get_NY_GRID();
    int dimz = get_NZ_GRID();
    size = dimx * dimy * dimz;
    int dist_size = get_PX0_GRID() *  get_PY0_GRID() * get_PZ0_GRID();

    for (isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        sp = &ct.sp[isp];
        std::complex<double> *saved_beta = (std::complex<double> *)sp->forward_beta;
        sp->forward_beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * dist_size * sp->num_projectors);
        std::complex<double> *fptr = (std::complex<double> *)sp->forward_beta;

        for(int iproj = 0;iproj < sp->num_projectors;iproj++) 
        {
            for(int i = 0;i < dimx;i++) {
                for(int j = 0;j < dimy;j++) {
                    for(int k = 0;k < dimz;k++) {
                        bool map = (i >= ilo) && (i < ihi) && (j >= jlo) && (j < jhi) && (k >= klo) && (k < khi);
                        if(map) {
                           int idx1 = (i - ilo) * (jhi - jlo) * (khi - klo) + (j - jlo) * (khi - klo) + (k - klo);
                           int idx2 = i * dimy * dimz + j * dimz + k;
                           fptr[iproj * dist_size + idx1] = saved_beta[iproj * size + idx2]; 
                        }
                    }
                }
            }
        }

    }

}                               /* end InitWeight */
