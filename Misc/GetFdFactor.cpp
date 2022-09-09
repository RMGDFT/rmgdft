
#include <float.h>
#include <math.h>
#include "main.h"
#include "Atomic.h"
#include "Pw.h"
#include "transition.h"


void GetFdFactor(void)
{
    FiniteDiff FD(&Rmg_L);
    std::complex<double> I_t(0.0, 1.0);

    int nlxdim = get_NX_GRID();
    int nlydim = get_NY_GRID();
    int nlzdim = get_NZ_GRID();
    int pbasis = Rmg_G->get_P0_BASIS(1);
    int norbitals = 0;

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
        sp.fd_factors.clear();
        sp.fd_fke.clear();

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
            norbitals++;

            // Get the FFT laplacian and compute the kinetic energy as our gold standard
            FftLaplacianCoarse(orbital, work);

            double fft_ke = 0.0;
            for(int idx=0;idx<pbasis;idx++) fft_ke += orbital[idx] * work[idx];
            fft_ke = -0.5*fft_ke*get_vel();
            MPI_Allreduce(MPI_IN_PLACE, &fft_ke, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

            double c2 = 0.0;
            cvals.clear();
            diffs.clear();
            // Just doing a linear fit with 2 points for now
            for(int j=0;j < 2;j++)
            {
                FD.cfac[0] = c2;
                if(ct.kohn_sham_fd_order == 8) ApplyAOperator (orbital, work, kvec);
                if(ct.kohn_sham_fd_order == 10) ApplyLaplacian (orbital, work, 10, "Coarse");
                double fd_ke = 0.0;
                for(int idx=0;idx<pbasis;idx++) fd_ke += orbital[idx] * work[idx];
                fd_ke = -0.5*fd_ke*get_vel();
                MPI_Allreduce(MPI_IN_PLACE, &fd_ke, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
                if(pct.gridpe == 0) printf("LLLL  %e   %e\n",c2, fft_ke - fd_ke);
                cvals.push_back(c2);
                diffs.push_back(fft_ke - fd_ke);
                c2 += 1.0; 
            }

            double m = (diffs[1] - diffs[0])/(cvals[1] - cvals[0]);
            double x_int = - diffs[0] / m;
            //if(x_int < 0.0) x_int = 0.0;
            //if(x_int > 1.0) x_int = 1.0;

            if(ct.verbose && pct.gridpe==0)printf("IP=%d M = %e  %e  %e\n",ip,m,x_int,diffs[0]);
            sp.fd_factors.push_back(x_int);
            // if the difference between fft and fd is very small then it likely
            // means the orbital is soft compared to the grid and we are just
            // fitting noise so don't include it in the computations since it
            // will have no impact on the total energy.
            if(fabs(m) > 1.0e-5)
                sp.fd_fke.push_back(fft_ke);
            else
                sp.fd_fke.push_back(0.0); // will cause it to not be included
        }
    } 

    // Loop over ions
    double newcfac = 0.0, kesum=0.0;
    for(auto& Atom : Atoms)
    {
        for(size_t i=0;i < Atom.Type->fd_factors.size();i++)
        {
            // Weight by the kinetic energy of the orbitals and the occupations
            newcfac += Atom.Type->fd_factors[i] * Atom.Type->fd_fke[i] * occ_weight[i];
            kesum += Atom.Type->fd_fke[i] * occ_weight[i];
        }
    }

    // If extrememly well converged then nothing to do here
    if(kesum == 0.0)
        newcfac = 0.0;
    else
        newcfac /= kesum;

    if(ct.verbose && pct.gridpe == 0) printf("NEWCFAC = %f\n",newcfac);
    FD.cfac[0] = newcfac;

    delete [] work;
    fftw_free (gbptr);
    fftw_free (beptr);
    delete [] orbital;
    delete [] fftw_phase;
}
