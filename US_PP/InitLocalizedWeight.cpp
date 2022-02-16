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


/*This sets loop over species does forward fourier transform, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void SPECIES::InitLocalizedWeight (void)
{

    fftw_complex *in, *out, *betaptr;
    std::complex<double> *phaseptr;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;

    
    RmgTimer RT0("Weight");
    // get tot number of projectors and their information

    int nldim_max = ct.max_nldim;
    
    RmgTimer *RT1= new RmgTimer("Weight: phase and set");

    int size = this->nldim * this->nldim * this->nldim;
    this->phase = new fftw_complex[size * ct.num_kpts_pe];
    phaseptr = (std::complex<double> *)this->phase;
    GetPhaseSpecies(this, phaseptr);
    /*Loop over all betas to calculate num of projectors for given species */
    int prjcount = 0;
    for (int ip = 0; ip < this->nbeta; ip++)
    {

        for(int m = 0; m < 2*this->llbeta[ip]+1; m++)
        {
            proj.ip = ip;
            proj.l = this->llbeta[ip];
            proj.m = m;
            proj.proj_index = prjcount;
            proj_iter.push_back(proj);
            prjcount++;
        }
    }

    size = this->nldim * this->nldim * this->nldim;
    this->num_projectors = prjcount;

    /*This array will store forward fourier transform on the coarse grid for all betas */
    if(this->forward_beta) fftw_free(this->forward_beta);
    this->forward_beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_projectors * size * ct.num_kpts_pe);

    if (this->forward_beta == NULL)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";


    delete RT1;

    int xdim = std::max(nldim_max, get_NX_GRID() );
    int ydim = std::max(nldim_max, get_NY_GRID() );
    int zdim = std::max(nldim_max, get_NZ_GRID() );
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * xdim * ydim * zdim);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * xdim * ydim * zdim);

    if(!in || !out)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";


    RmgTimer *RT3= new RmgTimer("Weight: proj cal");
    for(int iproj = pct.gridpe; iproj < this->num_projectors; iproj+=pct.grid_npes)
    {
        proj = proj_iter[iproj];

        // if this->nldim > get_NX_GRID, folding of neighbor cells are needed. 
        xdim = this->nldim;
        ydim = this->nldim;
        zdim = this->nldim;
        size = xdim * ydim * zdim;

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            phaseptr = (std::complex<double> *) &this->phase[kpt * this->nldim * this->nldim * this->nldim];
            betaptr = &this->forward_beta[kpt *this->num_projectors *size + proj.proj_index * size];
            InitWeightOne(this, betaptr, phaseptr, proj.ip, proj.l, proj.m);
        }

    }                           /* end for */

    fftw_free(out);
    fftw_free(in);
    delete RT3;
    RmgTimer *RT4= new RmgTimer("Weight: bcast");

    int root;
    for(int iproj = 0; iproj < this->num_projectors; iproj++)
    {
        proj = proj_iter[iproj];
        root = iproj % pct.grid_npes;
        size = this->nldim * this->nldim * this->nldim;

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            betaptr = &this->forward_beta[kpt *this->num_projectors *size + proj.proj_index * size];
            MPI_Bcast(betaptr, 2*size, MPI_DOUBLE, root, pct.grid_comm);
        }
    }

    delete [] this->phase;


    delete RT4;

}
