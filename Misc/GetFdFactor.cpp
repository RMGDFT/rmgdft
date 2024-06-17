
#include <float.h>
#include <math.h>
#include <boost/math/tools/minima.hpp>
#include <functional>
#include <boost/bind/bind.hpp> 
#include "main.h"
#include "Atomic.h"
#include "Pw.h"
#include "Lattice.h"
#include "transition.h"
#include "GlobalSums.h"

/*
   This routine is used to compute adaptive factors for the
   kinetic energy and the charge density. The procedure used
   for the kinetic energy is based on the method described in

   Adaptive finite differencing in high accuracy electronic structure calculations
   E. L. Briggs, Wenchang Lu & J. Bernholc 
   npj Computational Materials volume 10, Article number: 17 (2024) 

   but is modified to use the brent algorithm to find the value of M
   which minimizes the error in the kinetic energy. 


*/

using boost::math::tools::brent_find_minima;


// Evaluates a polynomial of order up to 8 using coefficients in coeffs
double eval_poly(double x, int order, std::array<double, 8> &coeffs)
{
    double fv = 0.0;
    for(int ik = order;ik > 0;ik--)
    {   
        fv = x*(fv + coeffs[ik]);
    }
    fv += coeffs[0];
    return fv;
}

// Evaluates the total kinetic energy error for use in the brent find minima routine
//
// x is the mixing M from the NPJ paper.
//
double eval_ke_error(double x)
{
    double kerr = 0.0;
    for(auto& Atom : Atoms)
    {
        for(size_t i=0;i < Atom.Type->fd_coeffs.size();i++)
        {
            kerr += Atom.Type->occ_weight[i] * eval_poly(x, 2, Atom.Type->fd_coeffs[i]);
        }
    }
    return kerr;
}


// Evaluates the total rho error for use in the brent find minima routine
//
// x is the mixing parameter for the higher and lower order operators
//
double eval_rho_error(double x)
{
    double err = 0.0;
    for(auto& Atom : Atoms)
    {
        for(size_t i=0;i < Atom.Type->fd_coeffs.size();i++)
        {
            err += Atom.Type->occ_weight[i] * eval_poly(x, 4, Atom.Type->pd_coeffs[i]);
        }
    }
    return err;
}

void SetCfacs(double *cf, double val)
{
    for(int i=0;i < 13;i++) cf[i] = val;
}

// Computes the kinetic energy of an orbital on the wavefunction grid
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

// Computes a measure of the effect on the charge density from
// interpolation errors. 
//
// x1 = FFT interpolated orbital on fine grid
// x2 = Prolong interpolated orbital on fine grid
//
double ComputeRhoGoodness(double *x1, double *x2, int pbasis)
{
    double rg1 = 0.0;
    for(int idx=0;idx<pbasis;idx++)
    {
        double t1 = (x1[idx] - x2[idx]);
        rg1 += std::pow(t1*t1, 1.0/3.0);
    }
    rg1 *= get_vel_f();
    MPI_Allreduce(MPI_IN_PLACE, &rg1, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    return rg1;
}


// Main routine which is only for gamma right now but may be extended
// to compute a different value for different k-points at some time
// in the future.
void GetFdFactor(int kpt)
{
    namespace ph = std::placeholders;

    // For the boost brent algorithm
    const int double_bits = std::numeric_limits<double>::digits;
    FiniteDiff FD(&Rmg_L);

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
    std::vector<double> cvals;
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    double *work = new double[pbasis];
    double *pwork1 = new double[fpbasis];
    double *pwork2 = new double[fpbasis];
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
        sp.pd_mins.clear();
        sp.occ_weight.clear();
        sp.fd_coeffs.clear();
        sp.pd_coeffs.clear();

        // Set up an occupation weight array
        for (int ip = 0; ip < sp.num_atomic_waves; ip++)
        {
            // This is here since we have forward beta only for occupied orbitals.
            // If that changes then this must change.
            if(sp.atomic_wave_oc[ip] > 0.0)
             {
                for(int m = 0; m < 2*sp.atomic_wave_l[ip]+1; m++)
                {
                    sp.occ_weight.push_back(sp.atomic_wave_oc[ip] / (2*sp.atomic_wave_l[ip]+1));
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
            std::vector<double> yarr1, yarr2;
            int N = 10;
            double dx = 4.0 / (double)N;
            // In this loop we generate data for the polynomial fits to use
            for(int j=0;j < N;j++)
            {

                SetCfacs(FD.cfac, c2);
                ApplyAOperator (orbital, work, kvec);
                double fd_ke = ComputeKineticEnergy(orbital, work, pbasis);
                if(ct.verbose && pct.gridpe == 0) 
                    fprintf(ct.logfile, "FFT-FD  %e   %e   %e   %e\n",c2, fft_ke, fd_ke, fft_ke - fd_ke);
                cvals.push_back(c2);
                yarr1.push_back((fft_ke - fd_ke)*(fft_ke - fd_ke));

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

            // Quadratic fit for KE.
            int order = 2;
            MPI_Barrier(MPI_COMM_WORLD);
            std::array<double, 8> coeffs;

            SimplePolyFit(cvals.data(), yarr1.data(), N, order, coeffs.data());
            std::pair<double, double> kmin;
            kmin = brent_find_minima(std::bind(eval_poly, ph::_1, order, coeffs), 0.0, 4.0, double_bits);
            if(ct.verbose && pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                printf("ke params  %14.8e  %14.8e\n", kmin.first, kmin.second);
            sp.fd_mins.push_back(kmin.first);
            sp.fd_coeffs.push_back(coeffs);

            // Quartic fit for rho
            MPI_Barrier(MPI_COMM_WORLD);
            order = 4;
            SimplePolyFit(cvals.data(), yarr2.data(), N, order, coeffs.data());
            sp.pd_coeffs.push_back(coeffs);

            std::pair<double, double> rhomin;
            rhomin = brent_find_minima(std::bind(eval_poly, ph::_1, order, coeffs), 0.0, 4.0, double_bits);
            if(ct.verbose && pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                printf("rho params  %14.8e  %14.8e\n",rhomin.first, rhomin.second);
            sp.pd_mins.push_back(rhomin.first);
        }
    }

    // Now we compute the M value that minimizes the total error in the KE.
    // or if the user entered a non-zero value in the input file we use that
    if(ct.afd_cfac > 0.0)
    {
        FD.cfac[0] = ct.afd_cfac;
        FD.cfac[1] = ct.afd_cfac;
    }
    else
    {
        std::pair<double, double> mmin;
        mmin = brent_find_minima(eval_ke_error, 0.0, 4.0, double_bits);
        FD.cfac[0] = mmin.first;
        FD.cfac[1] = mmin.first;
    }

    if(ct.verbose && pct.gridpe == 0 && pct.spinpe == 0 && kpt == 0)
    {
        fprintf(ct.logfile,"NEWCFAC = %f  %f\n",FD.cfac[0], FD.cfac[1]);
    }

    if(FD.cfac[0] < 0.0)
    {
        rmg_error_handler (__FILE__, __LINE__, 
            "CFAC < 0.0. This probably indicates an error in the cell setup:\n");
    }

    std::pair<double, double> cmin;
    cmin = brent_find_minima(eval_rho_error, 0.0, 4.0, double_bits);

    if(ct.use_cmix && ct.prolong_order)
    {
        if(ct.norm_conserving_pp)
        {
            ct.cmix = cmin.first;
        }
        else
        {
            if(pct.gridpe==0 && pct.spinpe==0 && kpt==0)
                printf("Notice: Adaptive interpolation is disabled for USPP and NLCC pseudopotentials.\n");
        }
    }

    delete [] pwork2;
    delete [] pwork1;
    delete [] work;
    fftw_free (gbptr);
    fftw_free (beptr);
    delete [] orbital;
    delete [] fftw_phase;
}
