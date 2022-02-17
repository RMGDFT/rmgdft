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


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void SPECIES::InitLocalizedOrbital (void)
{
    fftw_complex *betaptr;
    std::complex<double> *phaseptr;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;

    
    RmgTimer RT0("Weight");
    // get tot number of orbitals and their information

    int tot_orbitals = 0;
    
    RmgTimer *RT1= new RmgTimer("Weight: phase and set");
    int orbital_count = 0;

    int size = this->nldim * this->nldim * this->nldim;
    this->phase = new fftw_complex[size * ct.num_kpts_pe];
    phaseptr = (std::complex<double> *)this->phase;
    GetPhaseSpecies(this, phaseptr);

    /*Loop over all atomic waves to calculate num for given species */
    orbital_count = 0;
    for (int ip = 0; ip < this->num_atomic_waves; ip++)
    {
        for(int m = 0; m < 2*this->atomic_wave_l[ip]+1; m++)
        {
            proj.ip = ip;
            proj.l = this->atomic_wave_l[ip];
            proj.m = m;
            proj.proj_index = orbital_count;
            proj_iter.push_back(proj);
            orbital_count++;
        }
    }
    this->num_orbitals = orbital_count;

    size = this->nldim * this->nldim * this->nldim;

    if(this->forward_orbital) fftw_free(this->forward_orbital);
    this->forward_orbital = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_orbitals * size * ct.num_kpts_pe);

    if (this->forward_orbital == NULL)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    tot_orbitals += orbital_count;


    delete RT1;

    int xdim = std::max(this->nldim, get_NX_GRID() );
    int ydim = std::max(this->nldim, get_NY_GRID() );
    int zdim = std::max(this->nldim, get_NZ_GRID() );
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * xdim * ydim * zdim);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * xdim * ydim * zdim);

    if(!in || !out)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    RmgTimer *RT3= new RmgTimer("Weight: proj cal");
    for(int iproj = pct.gridpe; iproj < tot_orbitals; iproj+=pct.grid_npes)
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
            betaptr = &this->forward_orbital[kpt *this->num_projectors *size + proj.proj_index * size];
            InitWeightOne(this, betaptr, phaseptr, proj.ip, proj.l, proj.m);
        }

    }                           /* end for */

    fftw_free(out);
    fftw_free(in);
    delete RT3;
    RmgTimer *RT4= new RmgTimer("Orbital: bcast");

    int root;
    for(int iproj = 0; iproj < tot_orbitals; iproj++)
    {
        proj = proj_iter[iproj];
        root = iproj % pct.grid_npes;
        size = this->nldim * this->nldim * this->nldim;

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            betaptr = &this->forward_orbital[kpt *this->num_projectors *size + proj.proj_index * size];
            MPI_Bcast(betaptr, 2*size, MPI_DOUBLE, root, pct.grid_comm);
        }
    }

    delete [] this->phase;

    delete RT4;

} // end InitOrbital
