#include <iostream>     // std::cout
#include <algorithm>    // std::rotate
#include <vector>       // std::vector
#include <math.h>       // std::vector
#include <mpi.h>       // std::vector
#include "blas.h"
#include "GlobalSums.h"
#include "transition.h"
#include "RmgParallelFft.h"
#include "RmgException.h"

void Precond_drho(double *drho)
{

    // Kerker preconditioning parameter q0 in unit of au^-1;
    double q0 = ct.drho_q0;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;

    double beta = ct.resta_beta;
    double ktf = std::pow(3.0 * PI * ct.nel / Rmg_L.get_omega(), 1.0/3.0);
    double kappa = sqrt(ktf / PI) / 2.0;
    double gamma = sinh(kappa*beta) / kappa / beta;
    size_t pbasis = fine_pwaves->pbasis;
    std::vector<std::complex<double>> c_fm;
    c_fm.resize(pbasis);

    //printf("GGGGGG  %14.8e  %14.8e  %14.8e  %14.8e\n",ktf, beta, kappa, gamma); 
    for(int is = 0; is < ct.noncoll_factor * ct.noncoll_factor; is++)
    {
        for(size_t ig=0; ig < pbasis; ig++) c_fm[ig] = std::complex<double>(drho[is*pbasis + ig], 0.0);
        fine_pwaves->FftForward(c_fm.data(), c_fm.data());

        if(ct.drho_precond_type == PRECOND_RESTA)
        {
            for(size_t ig=0;ig < pbasis;ig++) {
                double g2 = fine_pwaves->gmags[ig] * tpiba2;
                if( g2 < 1.0e-5) g2 =  tpiba2;

                double a1 = kappa*kappa*sin(sqrt(g2)*beta);
                double a2 = gamma*sqrt(g2)*beta;
                double alpha = a1/a2 + g2;
                alpha = alpha / (g2 + kappa*kappa);
                c_fm[ig] = c_fm[ig] * alpha;
            }
        }
        else
        {
            for(size_t ig=0;ig < pbasis;ig++) {
                // drho is the charge density difference, sums of them always equal to zero
                // so there is no G= 0 term.
                // it is not true for spin polarized case or rho_x(yz) in SOC 
                // so we choose the smallest g2 for G=0 term
                double g2 = fine_pwaves->gmags[ig] * tpiba2;
                if( g2 < 1.0e-5) g2 =  tpiba2;

                double alpha = g2/(g2+ q0 * q0);
                c_fm[ig] = c_fm[ig] * alpha;
                // if(!fine_pwaves->gmask[ig]) c_fm[ig]  = 0.0;
                //c_fm[ig] = c_fm[ig] * g2/(g2+this->ktf * this->ktf);
            }
        }

        fine_pwaves->FftInverse(c_fm.data(), c_fm.data());
        for(size_t i = 0;i < pbasis;i++) drho[is*pbasis + i] = std::real(c_fm[i])/(double)fine_pwaves->global_basis;
    }
}
