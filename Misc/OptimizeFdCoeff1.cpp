
#include <float.h>
#include <math.h>
#include "main.h"
#include "Atomic.h"
#include "Pw.h"
#include "Lattice.h"
#include "transition.h"
#include "llbfgs.hpp"
#include "FDOpt.h"

void compute_der(int num_orb, int dimx, int dimy, int dimz, int pbasis, int sbasis, 
        int num_coeff, int order, double *orbitals, double * orbitals_b, double *psi_psin);
void compute_coeff_grad(int num_orb, std::vector<double> &ke_fft, std::vector<double> &ke_fd, 
        double *psi_psin, std::vector<double> &occ_weight, std::vector<double>& coeff_grad);
double ComputeKineticEnergy(double *x, double *lapx, int pbasis);

int FDOpt::pbasis;
int FDOpt::sbasis;
int FDOpt::order;
int FDOpt::num_orb;
int FDOpt::num_coeff;
std::vector<double> FDOpt::ke_fft;
std::vector<double> FDOpt::ke_fd;
std::vector<double> FDOpt::occ_weight;
std::vector<double> FDOpt::coeff;
std::vector<double> FDOpt::coeff_grad;
double *FDOpt::orbitals;
double *FDOpt::orbitals_b;
double *FDOpt::psi_psin;
double *FDOpt::work;
double FDOpt::kvec[3];


FDOpt::FDOpt(void)
{
    FiniteDiff FD(&Rmg_L);
    int nlxdim = get_NX_GRID();
    int nlydim = get_NY_GRID();
    int nlzdim = get_NZ_GRID();
    int dimx  = Rmg_G->get_PX0_GRID(1);
    int dimy  = Rmg_G->get_PY0_GRID(1);
    int dimz  = Rmg_G->get_PZ0_GRID(1);
    kvec[0] = 0.0;
    kvec[1] = 0.0;
    kvec[2] = 0.0;
    order = LC->Lorder;
    pbasis = Rmg_G->get_P0_BASIS(1);
    sbasis = (dimx + order) * (dimy + order) * (dimz + order);
    work = new double[sbasis];
 
    //  determine total number of orbital to be included, sum of atomic orbitals for each speices
    num_orb = 0;
    for (auto& sp : Species)
    {
        num_orb += sp.num_orbitals;
        sp.num_atoms = 0;
    }
    for(auto& Atom : Atoms)
    {
        Atom.Type->num_atoms +=1;
    }

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];

    // determine and initialize coeffients
    double c2 = FD.cfac[0];
    double c1 = 1.0+c2;
//  LC_6 is 2 orders lower than LC, not necessary to be 6 order
    ct.alt_laplacian = 0;
    for (int ax = 0; ax < 13; ax++)
    {
        if(!LC->include_axis[ax]) continue;
        LC->plane_centers[ax] = c1 * LC->plane_centers[ax] - c2 * LC_6->plane_centers[ax];
        coeff.push_back(LC->axis_lc[ax][0]);
        for(int i = 1; i < LC->Lorder/2; i++)
        {
            coeff.push_back(c1*LC->axis_lc[ax][i] - c2 * LC_6->axis_lc[ax][i-1]);
        } 
    }
    coeff_grad.resize(coeff.size());
    num_coeff = coeff.size();

    orbitals = new double[num_orb * pbasis];
    orbitals_b = new double[num_orb * sbasis];
    psi_psin = new double[num_orb * coeff.size()];
    ke_fft.resize(num_orb);
    ke_fd.resize(num_orb);

    double tot_occ = 0.0;
    for (auto& sp : Species)
    {
        // Set up an occupation weight array
        for (int ip = 0; ip < sp.num_atomic_waves; ip++)
        {
            // This is here since we have forward beta only for occupied orbitals.
            // If that changes then this must change.
            if(sp.atomic_wave_oc[ip] > 0.0)
            {
                for(int m = 0; m < 2*sp.atomic_wave_l[ip]+1; m++)
                {
                    occ_weight.push_back(sp.atomic_wave_oc[ip] / (2*sp.atomic_wave_l[ip]+1) * sp.num_atoms);
                    tot_occ += sp.atomic_wave_oc[ip] / (2*sp.atomic_wave_l[ip]+1) * sp.num_atoms;
                }
            }
        }
    }

    for (auto &occ_w : occ_weight) occ_w /= tot_occ;
    if(occ_weight.size() - num_orb != 0)
    {
        printf("\n occ_weigh size %d != num_orb %d", (int)occ_weight.size(), num_orb);
    }

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
    int orb_idx = 0;
    int images = LC->Lorder/2;
    for (auto& sp : Species)
    {

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

            for (int idx = 0; idx < pbasis; idx++) orbitals[orb_idx*pbasis + idx] = std::real(beptr[idx]);

            Rmg_T->trade_imagesx (&orbitals[orb_idx*pbasis], &orbitals_b[orb_idx*sbasis], dimx, dimy, dimz, images, FULL_TRADE);

            /*Advance the fortward transform pointers */
            fptr += pbasis;

            // Get the FFT laplacian and compute the kinetic energy as our gold standard
            FftLaplacianCoarse(&orbitals[orb_idx*pbasis], work);
            ke_fft[orb_idx] = ComputeKineticEnergy(&orbitals[orb_idx*pbasis], work, pbasis);
            orb_idx++;
        }
    }

    fftw_free (gbptr);
    fftw_free (beptr);
    delete [] fftw_phase;

    // calculte <psi(ijk)| psi(ijk + neighbor)> which will never change with different coeffients.
    compute_der(num_orb, dimx, dimy, dimz, pbasis, sbasis, num_coeff, order, orbitals, orbitals_b, psi_psin);
} // end FDOpt::FDOpt

FDOpt::~FDOpt(void)
{
    delete [] work;
}


double FDOpt::evaluate(void *, 
                       const double *x,
                       double *g,
                       const int n)
{
        double ke_diff2=1.0;
        int icoeff = 0;
        for (int ax = 0; ax < 13; ax++)
        {
            if(!LC->include_axis[ax]) continue;
            LC->plane_centers[ax] = 0.0;
            for(int i = 0; i < LC->Lorder/2; i++)
            {
                LC->axis_lc[ax][i] = coeff[icoeff];
                LC->plane_centers[ax] += -coeff[icoeff] * 2.0;
                icoeff++;
            } 
        }
        ke_diff2 = 0.0;
        for(int iorb = 0; iorb < num_orb; iorb++)
        {
            ApplyAOperator (&orbitals[iorb*pbasis], work, kvec);
            ke_fd[iorb] = ComputeKineticEnergy(&orbitals[iorb*pbasis], work, pbasis);
            ke_diff2 += (ke_fd[iorb] - ke_fft[iorb]) *(ke_fd[iorb] - ke_fft[iorb]) * occ_weight[iorb];
        }
        compute_coeff_grad(num_orb, ke_fft, ke_fd, psi_psin, occ_weight, coeff_grad);
        for(int i=0;i < n;i++) g[i] = -coeff_grad[i];
        printf("ke_diff2 = %14.8e\n",ke_diff2);
        return ke_diff2;

} // end FDOpt::evaluate


double FDOpt::Optimize(void)
{
    double ke_diff2;
    llbfgs::lbfgs_parameter_t lbfgs_params;
    llbfgs::lbfgs_load_default_parameters(&lbfgs_params);
    int ret = llbfgs::lbfgs_optimize(num_coeff, coeff.data(), &ke_diff2, evaluate, NULL, NULL, this, &lbfgs_params);
    int icoeff = 0;
    for (int ax = 0; ax < 13; ax++)
    {
        if(!LC->include_axis[ax]) continue;
        LC->plane_centers[ax] = 0.0;
        for(int i = 0; i < LC->Lorder/2; i++)
        {
            LC->axis_lc[ax][i] = coeff[icoeff];
            LC->plane_centers[ax] += -coeff[icoeff] * 2.0;
            icoeff++;
        } 
        }
    return ke_diff2;
}