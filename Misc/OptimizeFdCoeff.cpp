
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
        double *psi_psin, std::vector<double> &occ_weight, std::vector<double>& coeff_grad, std::vector<double>& coeff);
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
//  LC_6 is 2 orders lower than LC, not necessary to be 6 order
    ct.alt_laplacian = false;
    for (int ax = 0; ax < 13; ax++)
    {
        double c2 = FD.cfac[ax];
        double c1 = 1.0+c2;
        if(!LC->include_axis[ax]) continue;
        LC->plane_centers[ax] = c1 * LC->plane_centers[ax] - c2 * LC_6->plane_centers[ax];
        LC->axis_lc[ax][0] *= c1;
        coeff.push_back(LC->axis_lc[ax][0]);
        for(int i = 1; i < LC->Lorder/2; i++)
        {
            LC->axis_lc[ax][i] = c1*LC->axis_lc[ax][i] - c2 * LC_6->axis_lc[ax][i-1];
            coeff.push_back(LC->axis_lc[ax][i]);
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
    FDOpt::ke_fft.clear();
    FDOpt::ke_fd.clear();
    FDOpt::occ_weight.clear();
    FDOpt::coeff.clear();
    FDOpt::coeff_grad.clear();
    delete [] FDOpt::psi_psin;
    delete [] FDOpt::orbitals_b;
    delete [] FDOpt::orbitals;
    delete [] work;
}

// Fixed for now
double FDOpt::stepbound(void *,
                         const double *xp,
                         const double *d,
                         const int n)
{
  return 0.5;
}

int FDOpt::progress(void *instance,
                         const double *x,
                         const double *g,
                         const double fx,
                         const double xnorm,
                         const double gnorm,
                         const double step,
                         int n,
                         int k,
                         int ls)
{
    if(ct.verbose && pct.gridpe==0)
    {
        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
    }
    return 0;
}


double FDOpt::evaluate(void *, 
                       const double *x,
                       double *g,
                       const int n)
{
        double ke_diff2=0.0, poly_diff2 = 0.0;
        int icoeff = 0;
        for (int ax = 0; ax < 13; ax++)
        {
            if(!LC->include_axis[ax]) continue;
            LC->plane_centers[ax] = 0.0;
            double dist = 1.0 / LC->plane_dists[ax];
            double x2sum = 0.0;
            double x4sum = 0.0;
            double x6sum = 0.0;
            for(int i = 0; i < LC->Lorder/2; i++)
            {
                double h2 = dist * (double)(LC->Lorder/2 - i);
                h2 = h2 * h2;
                x2sum += coeff[icoeff]*h2;
                x4sum += coeff[icoeff]*h2*h2;
                x6sum += coeff[icoeff]*h2*h2*h2;
                LC->axis_lc[ax][i] = coeff[icoeff];
                if(ct.verbose) printf("For axis = %d  Coeff = %20.12f\n",ax,coeff[icoeff]);
                LC->plane_centers[ax] += -coeff[icoeff] * 2.0;
                icoeff++;
            } 
            poly_diff2 += (double)num_orb*(x2sum - 1.0)*(x2sum - 1.0);
            poly_diff2 += (double)num_orb*x4sum*x4sum;
            poly_diff2 += (double)num_orb*x6sum*x6sum;
        }
        for(int iorb = 0; iorb < num_orb; iorb++)
        {
            ApplyAOperator (&orbitals[iorb*pbasis], work, kvec);
            ke_fd[iorb] = ComputeKineticEnergy(&orbitals[iorb*pbasis], work, pbasis);
            ke_diff2 += (ke_fd[iorb] - ke_fft[iorb]) *(ke_fd[iorb] - ke_fft[iorb]) * occ_weight[iorb];
        }
        compute_coeff_grad(num_orb, ke_fft, ke_fd, psi_psin, occ_weight, coeff_grad, coeff);
        for(int i=0;i < n;i++) g[i] = -coeff_grad[i];
        printf("ke_diff2 = %14.8e   %14.8e\n",ke_diff2, poly_diff2);
        return ke_diff2 + poly_diff2;

} // end FDOpt::evaluate


double FDOpt::Optimize(void)
{
    double ke_diff2=0.0;
    llbfgs::lbfgs_parameter_t lbfgs_params;
    llbfgs::lbfgs_load_default_parameters(&lbfgs_params);
    lbfgs_params.delta = 0.0;
    lbfgs_params.past = 5;
    lbfgs_params.line_search_type = 0;
    int ret = llbfgs::lbfgs_optimize(num_coeff, coeff.data(), &ke_diff2, evaluate, stepbound, progress, this, &lbfgs_params);

    if(ct.verbose && pct.gridpe==0)
        printf("FD Opt lbfgs return msg = %s\n", llbfgs::lbfgs_strerror(ret));

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



void compute_coeff_grad(int num_orb, std::vector<double> &ke_fft, std::vector<double> &ke_fd, double *psi_psin, 
        std::vector<double> &occ_weight, std::vector<double> &coeff_grad, std::vector<double>& coeff)
{
    std::fill(coeff_grad.begin(), coeff_grad.end(), 0.0);
    int num_coeff = coeff_grad.size();
    std::vector<double> tcoeff_grad;
    tcoeff_grad.resize(num_coeff);
    std::fill(tcoeff_grad.begin(), tcoeff_grad.end(), 0.0);
    double gnorm1 = 0.0;
    for(int iorb = 0; iorb < num_orb; iorb++)
    {
        for(int i = 0; i < num_coeff; i++)
        {
            coeff_grad[i] += 2.0 * (ke_fd[iorb] - ke_fft[iorb]) * psi_psin[iorb * num_coeff + i] * occ_weight[iorb] ;
            coeff_grad[i] += 2.0 * (ke_fd[iorb] - ke_fft[iorb]) * psi_psin[iorb * num_coeff + i] * occ_weight[iorb] ;
        }
    }
    for(int i=0;i < num_coeff;i++) gnorm1 += coeff_grad[i]*coeff_grad[i];

        int icoeff = 0;
        for (int ax = 0; ax < 13; ax++)
        {
            if(!LC->include_axis[ax]) continue;
            double dist = 1.0 / LC->plane_dists[ax];
            double x2sum = 0.0;
            double x4sum = 0.0;
            double x6sum = 0.0;
            for(int i = 0; i < LC->Lorder/2; i++)
            {
                double h2 = dist * (double)(LC->Lorder/2 - i);
                h2 = h2 * h2;
                x2sum += coeff[icoeff]*h2;
                x4sum += coeff[icoeff]*h2*h2;
                x6sum += coeff[icoeff]*h2*h2*h2;
                icoeff++;
            } 
            int icoeff1 = icoeff - LC->Lorder/2;
printf("PPPP  %e  %e  %e\n",x2sum,x4sum,x6sum);
            for(int i = 0; i < LC->Lorder/2; i++)
            {
                double h2 = dist * (double)(LC->Lorder/2 - i);
                h2 = h2 * h2;
                tcoeff_grad[icoeff1] -= 2.0 * (x2sum - 1.0) * h2;
                tcoeff_grad[icoeff1] -= 2.0 * (x4sum) * h2 * h2;
                tcoeff_grad[icoeff1] -= 2.0 * (x6sum) * h2 * h2 * h2;
                icoeff1++;
            }
        }

    double gnorm2 = 0.0;
    for(int i=0;i < num_coeff;i++) gnorm2 += coeff_grad[i]*coeff_grad[i];
//    gnorm2 = sqrt(gnorm2 / gnorm1);
    for(int i=0;i < num_coeff;i++) coeff_grad[i] += (double)num_orb*tcoeff_grad[i];

}

void compute_der(int num_orb, int dimx, int dimy, int dimz, int pbasis, int sbasis, int num_coeff, int order, double *orbitals, double * orbitals_b, double *psi_psin)
{
    for(int i = 0; i < num_orb * num_coeff; i++) psi_psin[i] = 0.0;

    int ixs = (dimy + order) * (dimz + order);
    int iys = (dimz + order);
    int izs = 1;
    
    for (int iorb = 0; iorb < num_orb; iorb++)
    {
        double *b = &orbitals[iorb * pbasis]; 
        double *a = &orbitals_b[iorb * sbasis]; 
        int icoeff = 0;

        // take care of 3 axis along lattice vectors.
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                double *A = &a[iy*iys + ix*ixs];
                double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                // z-direction is orthogonal to xy-plane and only requires increments/decrements along z
                // 0=x,1=y,2=z,3=xy,4=xz,5=yz,6=nxy,7=nxz,8=nyz
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    for(int inn = 1; inn < order/2+1; inn++)
                    {
                        psi_psin[iorb * num_coeff + inn-1]  +=  B[iz] * (A[iz + inn * ixs] +  A[iz - inn * ixs] - A[iz]);
                        psi_psin[iorb * num_coeff + order/2 + inn-1]  +=  B[iz] * (A[iz + inn * iys] +  A[iz - inn * iys] - A[iz]);
                        psi_psin[iorb * num_coeff + 2 * order/2 + inn-1 ]  +=  B[iz] * (A[iz + inn * izs] +  A[iz - inn * izs] - A[iz]);
                    }
                }

            }
        }


        // Add additional axes as required
        icoeff = 3 * order/2;
        if(LC->include_axis[3])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz + inn * ixs + inn * iys] +  
                                    A[iz - inn * ixs - inn * iys] - A[iz]);
                        }
                    }                   /* end for */
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[4])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz + inn * ixs + inn * izs] +  
                                    A[iz - inn * ixs - inn * izs] - A[iz]);
                        }                   /* end for */
                    }
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[5])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz + inn * iys + inn * izs] +  
                                    A[iz - inn * iys - inn * izs] - A[iz]);
                        }                   /* end for */
                    }
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[6])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz - inn * ixs + inn * iys] +  
                                    A[iz + inn * ixs - inn * iys] - A[iz]);
                        }                   /* end for */
                    }
                }
            }
            icoeff += order/2;

        }

        if(LC->include_axis[7])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz - inn * ixs + inn * izs] +  
                                    A[iz + inn * ixs - inn * izs] - A[iz]);
                        }                   /* end for */
                    }
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[8])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz - inn * iys + inn * izs] +  
                                    A[iz + inn * iys - inn * izs] - A[iz]);
                        }
                    }
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[9])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz - inn * ixs - inn * iys - inn * izs] +  
                                    A[iz + inn * ixs + inn * iys + inn * izs] - A[iz]);
                        }                   /* end for */
                    }
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[10])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz - inn * ixs - inn * iys + inn *izs] +  
                                    A[iz + inn * ixs + inn * iys - inn *izs] - A[iz]);
                        }                   /* end for */
                    }
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[11])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz + inn * ixs - inn * iys + inn *izs] +  
                                    A[iz - inn * ixs + inn * iys - inn *izs] - A[iz]);
                        }                   /* end for */
                    }
                }
            }
            icoeff += order/2;
        }

        if(LC->include_axis[12])
        {
            for (int ix = order/2; ix < dimx + order/2; ix++)
            {
                for (int iy = order/2; iy < dimy + order/2; iy++)
                {
                    double *A = &a[iy*iys + ix*ixs];
                    double *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                    for (int iz = order/2; iz < dimz + order/2; iz++)
                    {
                        for(int inn = 1; inn < order/2+1; inn++)
                        {
                            psi_psin[iorb * num_coeff + icoeff + inn-1]  +=  B[iz] * (A[iz + inn * ixs - inn * iys - inn *izs] +  
                                    A[iz - inn * ixs + inn * iys + inn *izs] - A[iz]);
                        }
                    }                   /* end for */
                }
            }
            icoeff += order/2;
        }

    }

    int idim = num_orb * num_coeff;
    MPI_Allreduce(MPI_IN_PLACE, psi_psin, idim, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    for(int idx = 0; idx < num_orb * num_coeff; idx++) psi_psin[idx] *= -0.5 * get_vel();
}
void FDOpt::Analyze_fft(int orb_index)
{
    double tpiba2 = 4.0 * PI * PI / (Rmg_L.celldm[0] * Rmg_L.celldm[0]);
    int nx = get_NX_GRID();
    int ny = get_NY_GRID();
    int nz = get_NZ_GRID();

    ct.alt_laplacian = false;
    double *one_orbital = &orbitals[orb_index * pbasis];
    double *work_glob = new double[nx * ny * nz];
    double *gmags_glob = new double[nx * ny * nz];
    DistributeToGlobal(coarse_pwaves->gmags, gmags_glob);
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    FftLaplacianCoarse(one_orbital, work);
    for(int idx = 0; idx < pbasis; idx++) gbptr[idx] = std::complex<double>(work[idx],0.0);
    coarse_pwaves->FftForward(gbptr, beptr);

    for(int idx = 0; idx < pbasis; idx++) work[idx] = std::real(beptr[idx]);
    DistributeToGlobal(work, work_glob);
    //  for xz plane
    if(pct.gridpe == 0)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            int idx = ix * ny * nz +ix ;
            printf("\n %10.3f %18.12e   FFT_xz", gmags_glob[idx] *tpiba2, work_glob[idx]);
        }
    }

    ApplyAOperator (one_orbital, work, kvec);
    for(int idx = 0; idx < pbasis; idx++) gbptr[idx] = std::complex<double>(work[idx],0.0);
    coarse_pwaves->FftForward(gbptr, beptr);

    for(int idx = 0; idx < pbasis; idx++) work[idx] = std::real(beptr[idx]);
    DistributeToGlobal(work, work_glob);
    //  for xz plane
    if(pct.gridpe == 0)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            int idx = ix * ny * nz +ix ;
            printf("\n %10.3f %18.12e    FD_xz", gmags_glob[idx] *tpiba2, work_glob[idx]);
        }
    }

    delete [] work_glob;
    fftw_free (gbptr);
    fftw_free (beptr);
}




