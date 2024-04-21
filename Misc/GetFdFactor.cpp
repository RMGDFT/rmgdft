
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

double ComputeRhoGoodness(double *x1, double *x2, int pbasis)
{
    double rg = 0.0;
    for(int idx=0;idx<pbasis;idx++)
    {
        rg += (x1[idx] - x2[idx])*(x1[idx] - x2[idx]);
    }
    rg *= get_vel_f();
    MPI_Allreduce(MPI_IN_PLACE, &rg, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    rg = sqrt(rg);
    return rg;
}


void GetFdFactor(int kpt)
{
    FiniteDiff FD(&Rmg_L);

    if(ct.afd_cfac > 0.0)
    {
        FD.cfac[0] = ct.afd_cfac;
        FD.cfac[1] = ct.afd_cfac;
        return;
    }

    std::complex<double> I_t(0.0, 1.0);

    int nlxdim = get_NX_GRID();
    int nlydim = get_NY_GRID();
    int nlzdim = get_NZ_GRID();
    int pbasis = Rmg_G->get_P0_BASIS(1);
    int ratio = Rmg_G->default_FG_RATIO;
    int fpbasis = Rmg_G->get_P0_BASIS(ratio);
    int pxdim = get_PX0_GRID();
    int pydim = get_PY0_GRID();
    int pzdim = get_PZ0_GRID();

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];
    double *orbital = new double[pbasis];
    std::vector<double> cvals, diffs;
    std::vector<double> occ_weight;
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    double *work = new double[pbasis];
    double *pwork1 = new double[fpbasis];
    double *pwork2 = new double[fpbasis];
    double *pwork3 = new double[fpbasis];
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
        sp.pd_mins.clear();

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
        std::complex<double> *fptr = (std::complex<double> *)sp.forward_orbital_gamma;

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
                    fprintf(ct.logfile, "FFT-FD  %e   %e   %e   %e\n",c2, fft_ke, fd_ke, fft_ke - fd_ke);
                cvals.push_back(c2);
                diffs.push_back(fft_ke - fd_ke);


                c2 += 1.0; 
            }
            // Setup data for adaptive finite differencing
            double m = (diffs[1] - diffs[0])/(cvals[1] - cvals[0]);
            double x_int = - diffs[0] / m;
            sp.fd_slopes.push_back(m);
            sp.fd_yint.push_back(diffs[0]);

            if(ct.verbose && pct.gridpe==0)fprintf(ct.logfile,"IP=%d M = %e  %e  %e\n",ip,m,x_int,diffs[0]);
            x_int = std::max(x_int, 0.0);
            sp.fd_xint.push_back(x_int);

            sp.fd_factor1.push_back(x_int);
            sp.fd_fke1.push_back(fft_ke);

            // Now we do adaptive interpolation
            // Get the FFT Prolong as our gold standard */
            if(ct.prolong_order > 2)
            {
                FftInterpolation(*Rmg_G, orbital, pwork1, ratio, false);
                // Do 10 and find minimum
                c2 = 0.0;
                double lastval = 0.0;
                for(int j=0;j <= 40;j++)
                {
                    Prolong P(ratio, ct.prolong_order, c2, *Rmg_T,  Rmg_L, *Rmg_G);
                    P.prolong(pwork2, orbital, ratio*pxdim, ratio*pydim, ratio*pzdim, 
                              pxdim, pydim, pzdim);

                    double rg_p = ComputeRhoGoodness(pwork1, pwork2, fpbasis);
                    if(j > 0 && rg_p >= lastval)
                    {
                        sp.pd_mins.push_back(c2);
                        break; 
                    }
                    lastval = rg_p;
                    if(ct.verbose && pct.gridpe == 0) 
                        fprintf(ct.logfile, "RHOG = %e  %e\n", c2, rg_p);
                    c2 += 0.05;
                }
            } 
        }
    }

    // Loop over ions
    double fweight=0.0, a=0.0;
    double pweight=0.0, p=0.0;
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

    ct.cmix = 1.0;
    if(0 && ct.prolong_order > 2)
    {
        for(auto& Atom : Atoms)
        {
            for(size_t i=0;i < Atom.Type->pd_mins.size();i++)
            {
                if(fabs(Atom.Type->pd_mins[i]) > 0.0 && fabs(Atom.Type->pd_mins[i]) < 2.0)
                {
                    p += Atom.Type->pd_mins[i] * occ_weight[i];
                    pweight += occ_weight[i];
                }
            }
        }
        if(pweight > 0)
        {
            if(ct.verbose && pct.gridpe == 0) 
                printf("cmix = %f\n", p/pweight);

            ct.cmix = p/pweight;
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

    if(ct.verbose && pct.gridpe == 0) fprintf(ct.logfile,"NEWCFAC = %f  %f\n",FD.cfac[0], FD.cfac[1]);

    if(FD.cfac[0] < 0.0)
    {
        rmg_error_handler (__FILE__, __LINE__, 
            "CFAC < 0.0. This probably indicates an error in the cell setup:\n");
    }

    delete [] pwork3;
    delete [] pwork2;
    delete [] pwork1;
    delete [] work;
    fftw_free (gbptr);
    fftw_free (beptr);
    delete [] orbital;
    delete [] fftw_phase;
}
