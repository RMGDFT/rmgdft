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


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void InitWeight (void)
{

    int ip, prjcount, isp, size, tot_proj;
    SPECIES *sp;
    fftw_plan p1;
    fftw_complex *in, *out;
    int maxl = 0;

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
        sp->num_projectors = prjcount;
        sp->forward_beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_projectors * size);
        if (sp->forward_beta == NULL)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

        tot_proj += prjcount;
    }


    for(int iproj = pct.gridpe; iproj < tot_proj; iproj+=pct.grid_npes)
    {
        proj = proj_iter[iproj];

        sp = &ct.sp[proj.species];
 
        size = sp->nldim * sp->nldim * sp->nldim;


        /*This is something we need to do only once per species, so do not use wisdom */
        in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
        out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);

        if(!in || !out)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

        p1 = fftw_plan_dft_3d (sp->nlfdim, sp->nlfdim, sp->nlfdim, in, out, FFTW_FORWARD, FFTW_MEASURE);

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

        MPI_Bcast(&sp->forward_beta[proj.proj_index * size], 2*size, MPI_DOUBLE, root, pct.grid_comm);
    }
        
}                               /* end get_weight */
