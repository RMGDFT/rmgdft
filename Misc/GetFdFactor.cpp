
#include <float.h>
#include <math.h>
#include "main.h"
#include "Atomic.h"
#include "Pw.h"
#include "Lattice.h"
#include "transition.h"


void SetCfacs(double *cf, double val)
{
    for(int i=0;i < 13;i++) cf[i] = val;
}

double ComputeKineticEnergy(double *x, double *lapx, int pbasis)
{
    double ke = 0.0;
    for(int idx=0;idx<pbasis;idx++)
    {
        ke += x[idx] * lapx[idx];
    }
    ke = -0.5*ke*get_vel();
    MPI_Allreduce(MPI_IN_PLACE, &ke, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    return ke;
}


void GetFdFactor(int kpt)
{
    FiniteDiff FD(&Rmg_L);
    std::complex<double> I_t(0.0, 1.0);

    int nlxdim = get_NX_GRID();
    int nlydim = get_NY_GRID();
    int nlzdim = get_NZ_GRID();
    int pbasis = Rmg_G->get_P0_BASIS(1);

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];
    double *orbital = new double[pbasis];
    std::vector<double> cvals, diffs;
    std::vector<double> occ_weight;
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    double *work = new double[pbasis];
    double vect[3], nlcrds[3], kvec[3];

    /* Find nlcdrs, vector that gives shift of ion from center of its ionic box */
    /* for delocalized case it's just half the cell dimensions */
    vect[0] = 0.5;
    vect[1] = 0.5;
    vect[2] = 0.5;
    kvec[0] = 0.0;
    kvec[1] = 0.0;
    kvec[2] = 0.0;

    // Loop over species
    for (auto& sp : Species)
    {
        sp.fd_factor1.clear();
        sp.fd_fke1.clear();
        sp.fd_slopes.clear();
        sp.fd_xint.clear();
        sp.fd_yint.clear();

        // Set up an occupation weight array
        for (int ip = 0; ip < sp.num_atomic_waves; ip++)
        {
            // This is here since we have forward beta only for occupied orbitals.
            // If that changes then this must change.
            if(sp.atomic_wave_oc[ip] > 0.0)
             {
                for(int m = 0; m < 2*sp.atomic_wave_l[ip]+1; m++)
                {
                    occ_weight.push_back(sp.atomic_wave_oc[ip] / (2*sp.atomic_wave_l[ip]+1));
                }
             }
        }
 
        /*The vector we are looking for should be */
        to_cartesian (vect, nlcrds);

        /*Calculate the phase factor */
        FindPhaseKpoint (kvec, nlxdim, nlydim, nlzdim, nlcrds, fftw_phase, false);

        /*Temporary pointer to the already calculated forward transform. */
        /* Need to fix up for kpoint parrallelization issues.  */
        std::complex<double> *fptr = (std::complex<double> *)sp.forward_orbital;
        fptr += kpt*sp.num_orbitals*pbasis;

        /* Loop over atomic orbitals */
        for (int ip = 0; ip < sp.num_orbitals; ip++)
        {
            /*Apply the phase factor */
            for (int idx = 0; idx < pbasis; idx++) gbptr[idx] = fptr[idx] * std::conj(fftw_phase[idx]);

            /*Do the backwards transform */
            coarse_pwaves->FftInverse(gbptr, beptr);

            for (int idx = 0; idx < pbasis; idx++) orbital[idx] = std::real(beptr[idx]);

            /*Advance the fortward transform pointers */
            fptr += pbasis;

            // Get the FFT laplacian and compute the kinetic energy as our gold standard
            FftLaplacianCoarse(orbital, work);
            double fft_ke = ComputeKineticEnergy(orbital, work, pbasis);
            double c2 = 0.0;
            cvals.clear();
            diffs.clear();
            // Linear fit with 2 points
            for(int j=0;j < 2;j++)
            {
                SetCfacs(FD.cfac, c2);
                ApplyAOperator (orbital, work, kvec);
                double fd_ke = ComputeKineticEnergy(orbital, work, pbasis);
                if(ct.verbose && pct.gridpe == 0) 
                    printf("FFT-FD  %e   %e   %e   %e\n",c2, fft_ke, fd_ke, fft_ke - fd_ke);
                cvals.push_back(c2);
                diffs.push_back(fft_ke - fd_ke);
                c2 += 1.0; 
            }
            double m = (diffs[1] - diffs[0])/(cvals[1] - cvals[0]);
            double x_int = - diffs[0] / m;
            sp.fd_slopes.push_back(m);
            sp.fd_yint.push_back(diffs[0]);

            if(ct.verbose && pct.gridpe==0)printf("IP=%d M = %e  %e  %e\n",ip,m,x_int,diffs[0]);
            x_int = std::max(x_int, 0.0);
            sp.fd_xint.push_back(x_int);

            sp.fd_factor1.push_back(x_int);
            sp.fd_fke1.push_back(fft_ke);
        }
    }

    // Loop over ions
    double fweight=0.0, a=0.0;
    for(auto& Atom : Atoms)
    {
        for(size_t i=0;i < Atom.Type->fd_slopes.size();i++)
        {
            if(fabs(Atom.Type->fd_slopes[i]) > 1.0e-8)
            {
                a += Atom.Type->fd_slopes[i] * occ_weight[i] * 
                     (Atom.Type->fd_xint[i] - Atom.Type->fd_yint[i]);
                fweight += Atom.Type->fd_slopes[i] * occ_weight[i];
            }
        }
    }

    // If extrememly well converged then nothing to do here
    if(fweight == 0.0)
    {
        FD.cfac[0] = 0.0;
        FD.cfac[1] = 0.0;
    }
    else
    {
        FD.cfac[0] = a / fweight;
        FD.cfac[1] = a / fweight;
    }

    if(ct.verbose && pct.gridpe == 0) printf("NEWCFAC = %f  %f\n",FD.cfac[0], FD.cfac[1]);

    if(FD.cfac[0] < 0.0)
    {
        rmg_error_handler (__FILE__, __LINE__, 
            "CFAC < 0.0. This probably indicates an error in the cell setup:\n");
    }

    delete [] work;
    fftw_free (gbptr);
    fftw_free (beptr);
    delete [] orbital;
    delete [] fftw_phase;
}
