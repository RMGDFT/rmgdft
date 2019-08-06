#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <complex>


#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgException.h"
#include "RmgShm.h"


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void InitOrbital (void)
{

    int ip, isp, size, tot_orbitals;
    SPECIES *sp;
    fftw_plan p2_forward, p2_backward;
    fftw_complex *in, *out, *betaptr;
    std::complex<double> *phaseptr;
    int xdim, ydim, zdim, nldim_max;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;

    
    RmgTimer RT0("Weight");
    // get tot number of projectors and their information

    nldim_max = 0;
    tot_orbitals = 0;
    
    RmgTimer *RT1= new RmgTimer("Weight: phase and set");
    for (isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        sp = &Species[isp];
        int orbital_count = 0;

        nldim_max = std::max(nldim_max, sp->nldim);
        size = sp->nldim * sp->nldim * sp->nldim;
        sp->phase = new fftw_complex[size * ct.num_kpts_pe];
        phaseptr = (std::complex<double> *)sp->phase;
        GetPhaseSpecies(sp, phaseptr);
        /*Loop over all betas to calculate num of projectors for given species */
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            for(int m = 0; m < 2*sp->llbeta[ip]+1; m++)
            {
                proj.species = isp;
                proj.ip = ip;
                proj.l = sp->atomic_wave_l[ip];
                proj.m = m;
                proj.proj_index = orbital_count;
                proj_iter.push_back(proj);
                orbital_count++;
            }
        }

        size = sp->nldim * sp->nldim * sp->nldim;

        
        orbital_count = 0;

        for (int ip = 0; ip < sp->num_atomic_waves; ip++) {
            if(sp->atomic_wave_oc[ip] > 0.0) {
                int l = sp->atomic_wave_l[ip];
                for (int m=0; m < 2*l+1; m++) {
                    orbital_count ++;
                }
            }
        }

        sp->num_orbitals = orbital_count;


        sp->forward_orbital = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_orbitals * size * ct.num_kpts_pe);

        if (sp->forward_orbital == NULL)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

        tot_orbitals += orbital_count;
    }

    delete RT1;

    xdim = std::max(nldim_max, get_NX_GRID() );
    ydim = std::max(nldim_max, get_NY_GRID() );
    zdim = std::max(nldim_max, get_NZ_GRID() );
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * xdim * ydim * zdim);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * xdim * ydim * zdim);

    if(!in || !out)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    RmgTimer *RT3= new RmgTimer("Weight: proj cal");
    for(int iproj = pct.gridpe; iproj < tot_orbitals; iproj+=pct.grid_npes)
    {
        proj = proj_iter[iproj];

        sp = &Species[proj.species];

        // if sp->nldim > get_NX_GRID, folding of neighbor cells are needed. 
        xdim = sp->nldim;
        ydim = sp->nldim;
        zdim = sp->nldim;
        size = xdim * ydim * zdim;

        p2_backward = fftw_plan_dft_3d (sp->nldim, sp->nldim, sp->nldim, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        p2_forward = fftw_plan_dft_3d (xdim, ydim, zdim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            phaseptr = (std::complex<double> *) &sp->phase[kpt * sp->nldim * sp->nldim * sp->nldim];
            betaptr = &sp->forward_orbital[kpt *sp->num_projectors *size + proj.proj_index * size];
            InitWeightOne(sp, betaptr, phaseptr, proj.ip, proj.l, proj.m, p2_forward, p2_backward);
        }

        fftw_destroy_plan (p2_forward);
        fftw_destroy_plan (p2_backward);

    }                           /* end for */

    fftw_free(out);
    fftw_free(in);
    delete RT3;
    RmgTimer *RT4= new RmgTimer("Orbital: bcast");

    int root;
    for(int iproj = 0; iproj < tot_orbitals; iproj++)
    {
        proj = proj_iter[iproj];
        sp = &Species[proj.species];
        root = iproj % pct.grid_npes;
        size = sp->nldim * sp->nldim * sp->nldim;

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            betaptr = &sp->forward_orbital[kpt *sp->num_projectors *size + proj.proj_index * size];
            MPI_Bcast(betaptr, 2*size, MPI_DOUBLE, root, pct.grid_comm);
        }
    }

    for (isp = 0; isp < ct.num_species; isp++)
    {
        sp = &Species[isp];
        delete sp->phase;
    }


    delete RT4;

} // end InitOrbital
