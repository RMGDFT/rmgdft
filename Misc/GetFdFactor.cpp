
#include <float.h>
#include <math.h>
#include <boost/math/tools/minima.hpp>
#include "main.h"
#include "Atomic.h"
#include "Pw.h"
#include "Lattice.h"
#include "transition.h"
#include "GlobalSums.h"

using boost::math::tools::brent_find_minima;



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
    double rg1 = 0.0;
    for(int idx=0;idx<pbasis;idx++)
    {
        double t1 = (x1[idx] - x2[idx]);
        rg1 += sqrt(t1*t1);
    }
    rg1 *= get_vel_f();
    MPI_Allreduce(MPI_IN_PLACE, &rg1, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    return rg1;
}


void GetFdFactor(int kpt)
{
    FiniteDiff FD(&Rmg_L);
    int is_core_correction = false;

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
        if(sp.nlccflag) is_core_correction = true;
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
 
            // Make sure the orbital is normalized to 1.0
            double snorm = 0.0;
            for(int idx=0;idx < pbasis;idx++) snorm += std::real(orbital[idx] * std::conj(orbital[idx]));
            GlobalSums(&snorm, 1, pct.grid_comm);
            snorm *= get_vel();
            snorm = 1.0 / sqrt(snorm);
            for(int idx=0;idx < pbasis;idx++) orbital[idx] *= snorm;
            if(ct.verbose && pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                printf("SNORM %d  %d  %14.8f\n",sp.num_orbitals,ip,snorm);

            /*Advance the fortward transform pointers */
            fptr += pbasis;

            // Get the FFT laplacian and compute the kinetic energy as our gold standard
            FftLaplacianCoarse(orbital, work);
            double fft_ke = ComputeKineticEnergy(orbital, work, pbasis);


            // Now do the FFT Interpolation to the fine grid to use as
            // a reference standard for optimizing the prolongation operator.
            FftInterpolation(*Rmg_G, orbital, pwork1, ratio, false);

            double c2 = 0.0;
            cvals.clear();
            diffs.clear();
            std::vector<double> yarr1, yarr2;
            int N = 10;
            double dx = 4.0 / (double)N;
            for(int j=0;j < N;j++)
            {

                SetCfacs(FD.cfac, c2);
                ApplyAOperator (orbital, work, kvec);
                double fd_ke = ComputeKineticEnergy(orbital, work, pbasis);
                if(ct.verbose && pct.gridpe == 0) 
                    fprintf(ct.logfile, "FFT-FD  %e   %e   %e   %e\n",c2, fft_ke, fd_ke, fft_ke - fd_ke);
                cvals.push_back(c2);
                diffs.push_back(fft_ke - fd_ke);
                yarr1.push_back(fft_ke - fd_ke);

                Prolong P(ratio, ct.prolong_order, c2, *Rmg_T,  Rmg_L, *Rmg_G);
                P.prolong(pwork2, orbital, ratio*pxdim, ratio*pydim, ratio*pzdim, 
                              pxdim, pydim, pzdim);

                double snorm_f = 0.0;
                for(int idx=0;idx < fpbasis;idx++) snorm_f += std::real(pwork2[idx] * std::conj(pwork2[idx]));
                GlobalSums(&snorm_f, 1, pct.grid_comm);
                snorm_f *= get_vel_f();
                snorm_f = 1.0 / sqrt(snorm_f);
                for(int idx=0;idx < fpbasis;idx++) pwork2[idx] *= snorm_f;

                double rhogood = ComputeRhoGoodness(pwork1, pwork2, fpbasis);
                if(ct.verbose && pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                    printf("rhogood  %14.8f   %14.8e\n",c2, rhogood);
                yarr2.push_back(rhogood);
                c2 += dx; 
            }
            if(ct.verbose && pct.gridpe==0 && pct.spinpe==0 && kpt==0)printf("rhogood  &&\n");

            // Setup data for adaptive finite differencing
            double m = (diffs[1] - diffs[0])/(cvals[1] - cvals[0]);
            double x_int = - diffs[0] / m;
            sp.fd_slopes.push_back(m);
            sp.fd_yint.push_back(diffs[0]);
            x_int = std::max(x_int, 0.0);
            sp.fd_xint.push_back(x_int);
            if(ct.verbose && pct.gridpe==0)
                fprintf(ct.logfile,"IP=%d M = %e  %e  %e\n",ip,m,x_int,diffs[0]);

            // Linear fit for KE
            MPI_Barrier(MPI_COMM_WORLD);
            static int korder = 1;
            static double coeffs[10];  // Fitted polynomial coefficients

            SimplePolyFit(cvals.data(), yarr1.data(), N, korder, coeffs);
            const int double_bits = std::numeric_limits<double>::digits;
            if(ct.verbose && pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                printf("ke params  %14.8f  %14.8f\n",coeffs[0], coeffs[1]);
            double fd_xint = -coeffs[0]/coeffs[1];  // Want y-intercept
            sp.fd_mins.push_back(fd_xint);

            // Quadratic fit for rho
            MPI_Barrier(MPI_COMM_WORLD);
            korder = 4;
            SimplePolyFit(cvals.data(), yarr2.data(), N, korder, coeffs);

            // Find minimum using brent algorithm from boost
            struct func_rho
            {
              double operator()(double const& x)
              {
                double fv = 0.0;
                for(int ik = korder;ik > 0;ik--)
                {
                    fv = x*(fv + coeffs[ik]);
                }
                fv += coeffs[0]; 
                return fv;
              }
            };
            std::pair<double, double> rhomin = brent_find_minima(func_rho(), 0.0, 4.0, double_bits);
            if(ct.verbose && pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                printf("rho params  %14.8e  %14.8e\n",rhomin.first, rhomin.second);
            sp.pd_mins.push_back(rhomin.first);
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
//                a += Atom.Type->fd_mins[i] * occ_weight[i];
//                fweight += occ_weight[i];
            }
            p += Atom.Type->pd_mins[i] * occ_weight[i];
            pweight += occ_weight[i];
        }
    }

    if(ct.use_cmix && ct.prolong_order > 2 && std::abs(pweight) > 1.0e-8)
    {
        if(is_core_correction == false && ct.norm_conserving_pp)
        {
            ct.cmix = p/pweight;
        }
        else
        {
            if(pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                printf("Notice: Adaptive interpolation is disabled for USPP and NLCC pseudopotentials.\n");
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
