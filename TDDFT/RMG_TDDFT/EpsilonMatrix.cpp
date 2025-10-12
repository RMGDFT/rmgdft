/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include <complex>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <libgen.h>

#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "blacs.h"
#include "RmgParallelFft.h"
#include "RmgException.h"


#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"



template void EpsilonMatrix<double>(Kpoint<double> **Kptr);
template void EpsilonMatrix<std::complex<double>>(Kpoint<std::complex<double>> **Kptr);

 template <typename KpointType>
void EpsilonMatrix (Kpoint<KpointType> **Kptr)
{

    // instead of calculating <psi|i nabla -k |psi>, we calculate <psi |nabla +i k|psi> 
    Kpoint<KpointType>  *kptr;
    kptr = Kptr[0];
    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    int num_states = kptr->nstates;
    int nb=std::min(num_states, ct.scalapack_block_factor);
    int pbasis = kptr->pbasis;
    int pbasis_noncol = pbasis * ct.noncoll_factor;
    int ix, iy, iz;
    Rmg_G->pe2xyz (pct.gridpe, &ix, &iy, &iz);
    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);

    int px0_grid = Rmg_G->get_PX0_GRID(1);
    int py0_grid = Rmg_G->get_PY0_GRID(1);
    int pz0_grid = Rmg_G->get_PZ0_GRID(1);
    double xoff = ix * px0_grid * hxgrid;
    double yoff = iy * py0_grid * hygrid;
    double zoff = iz * pz0_grid * hzgrid;

    double xtal[3], xcrt[3];

    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));
    //  alpha take care of i in moment operator
    KpointType alpha(vel);
    KpointType beta(0.0);

    
    std::vector<KpointType> Pmat[3];
    Pmat[0].resize(num_states * num_states,0.0);
    Pmat[1].resize(num_states * num_states,0.0);
    Pmat[2].resize(num_states * num_states,0.0);

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    trans_a = trans_c;

    //block_size = num_states;
    int num_blocks = (num_states + nb -1)/nb;
    // First time through allocate pinned memory for global_matrix1

    // 3 block matrix for px, py, pz operators
    KpointType *block_matrix;
    int retval1 = MPI_Alloc_mem(3*num_states * nb * sizeof(KpointType) , MPI_INFO_NULL, &block_matrix);
    if(retval1 != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
    }
    KpointType *block_matrix_x = block_matrix;
    KpointType *block_matrix_y = block_matrix_x + num_states * nb;
    KpointType *block_matrix_z = block_matrix_y + num_states * nb;

    double Emax = 50.0/Ha_eV;
    int Epoints = 5001;
    double delta_e = Emax/(Epoints-1);
    std::vector<double> epsilon[9];
    std::vector<double> epsilon_lorentzian[9];
    double eps = 1.0e-6;
    for(int i = 0; i < 9; i++)
    {
        epsilon[i].resize(Epoints, 0.0);
        epsilon_lorentzian[i].resize(Epoints, 0.0);
    }
    // ct.gaus_broad is in unit of eV
    double gaus_broad = ct.gaus_broad/Ha_eV;
    int negrid = 1+int(gaus_broad/delta_e *std::sqrt(-std::log(eps * gaus_broad*std::sqrt(PI))));
    double e2gaus = delta_e * delta_e /gaus_broad/gaus_broad;

    // delta function with n point away from center
    // std::exp(- n^2 * delt_e^2 /gaus_broad/gaus_broad) / gaus_broad /std::sqrt(PI);

    MPI_Barrier(MPI_COMM_WORLD);
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {  
        std::fill(Pmat[0].begin(),Pmat[0].end(), 0.0);
        std::fill(Pmat[1].begin(),Pmat[1].end(), 0.0);
        std::fill(Pmat[2].begin(),Pmat[2].end(), 0.0);
        kptr = Kptr[kpt];

        KpointType *ns = kptr->ns;
        KpointType *nv = kptr->nv;
        KpointType *newsint_local = kptr->newsint_local;



        // V|psi> is in tmp_arrayT
        KpointType *psi = kptr->orbital_storage;
        KpointType *psi_x = psi + num_states * pbasis_noncol;
        KpointType *psi_y = psi_x + nb * pbasis_noncol;
        KpointType *psi_z = psi_y + nb * pbasis_noncol;
        KpointType *psi_dev;
        if(kptr->psi_dev)
        {
            psi_dev = kptr->psi_dev;
        }
        else
        {
            psi_dev = psi;
        }


        for(int ib = 0; ib < num_blocks; ib++)
        {
            // upper triagle blocks only
            // ib : (nstates - ib * nb) * block_size matrix
            // ib = 0: nstates * block_size matrix
            // ib = 1: (nstates - nb) * block_size matrix
            // block index (ib, ib:num_blocks)
            // this_block_size will be nb for the first num_blocks-1 block, the last block could be smaller than nb


            int this_block_size, length_block;
            length_block = num_states - nb * ib;
            this_block_size = std::min(nb, length_block);

            int st_start = ib *nb;
            AppNls(kptr, newsint_local, kptr->Kstates[0].psi, nv, ns, st_start, this_block_size);

            for (int st1 = 0; st1 < this_block_size; st1++)
            {
                KpointType *psi1 = psi + (ib*nb + st1) * pbasis_noncol;
                KpointType *psi1_x = psi_x + st1 * pbasis_noncol;
                KpointType *psi1_y = psi_y + st1 * pbasis_noncol;
                KpointType *psi1_z = psi_z + st1 * pbasis_noncol;
                ApplyGradient(psi1, psi1_x, psi1_y, psi1_z, ct.force_grad_order, "Coarse");
                if(ct.noncoll)
                {
                    ApplyGradient(psi1+pbasis, psi1_x+pbasis, psi1_y+pbasis, psi1_z+pbasis, ct.force_grad_order, "Coarse");
                }

                if(!ct.is_gamma)
                {
                    std::complex<double> I_t(0.0, 1.0);
                    std::complex<double> *psi_C, *psi_xC, *psi_yC, *psi_zC;
                    psi_C = (std::complex<double> *) psi1;
                    psi_xC = (std::complex<double> *) psi1_x;
                    psi_yC = (std::complex<double> *) psi1_y;
                    psi_zC = (std::complex<double> *) psi1_z;
                    for(int i = 0; i < pbasis_noncol; i++)
                    {
                        psi_xC[i] += I_t *  kptr->kp.kvec[0] * psi_C[i];
                        psi_yC[i] += I_t *  kptr->kp.kvec[1] * psi_C[i];
                        psi_zC[i] += I_t *  kptr->kp.kvec[2] * psi_C[i];
                    }
                }

                for(int i = 0; i < px0_grid; i++)
                {
                    for(int j = 0; j < py0_grid; j++)
                    {
                        for(int k = 0; k < pz0_grid; k++)
                        {

                            xtal[0] = xoff + i * hxgrid;
                            xtal[1] = yoff + j * hygrid;
                            xtal[2] = zoff + k * hzgrid;
                            Rmg_L.to_cartesian(xtal, xcrt);

                            int idx = st1 * pbasis_noncol + i * py0_grid * pz0_grid + j * pz0_grid + k;
                            // nonlocal part need add the conj part
                            psi_x[idx] = 0.5 * psi_x[idx] + nv[idx] * xcrt[0] ;
                            psi_y[idx] = 0.5 * psi_y[idx] + nv[idx] * xcrt[1] ;
                            psi_z[idx] = 0.5 * psi_z[idx] + nv[idx] * xcrt[2] ;

                        } 
                    }
                }
            }

            RmgGemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_x, pbasis_noncol, psi_dev, 
                    pbasis_noncol, beta, block_matrix_x, this_block_size);
            BlockAllreduce((double *)block_matrix_x, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);
            RmgGemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_y, pbasis_noncol, psi_dev, 
                    pbasis_noncol, beta, block_matrix_y, this_block_size);
            BlockAllreduce((double *)block_matrix_y, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

            RmgGemm(trans_a, trans_n, this_block_size, num_states,  pbasis_noncol, alpha, psi_z, pbasis_noncol, psi_dev, 
                    pbasis_noncol, beta, block_matrix_z, this_block_size);
            BlockAllreduce((double *)block_matrix_z, (size_t)this_block_size * (size_t)num_states * (size_t)factor , pct.grid_comm);

            for(int j = 0; j < num_states; j++)
            {
                for(int i = 0; i < this_block_size; i++)
                {
                    Pmat[0][j * num_states + ib*nb + i] += block_matrix_x[j * this_block_size + i];
                    Pmat[1][j * num_states + ib*nb + i] += block_matrix_y[j * this_block_size + i];
                    Pmat[2][j * num_states + ib*nb + i] += block_matrix_z[j * this_block_size + i];
                    Pmat[0][(ib*nb + i) * num_states + j ] -= MyConj(block_matrix_x[j * this_block_size + i]);
                    Pmat[1][(ib*nb + i) * num_states + j ] -= MyConj(block_matrix_y[j * this_block_size + i]);
                    Pmat[2][(ib*nb + i) * num_states + j ] -= MyConj(block_matrix_z[j * this_block_size + i]);
                }
            }

        }

        int vmax_st = 0;
        int cmin_st = 0;

        double full_occ = 2.0;
        if(ct.nspin != 1) full_occ = 1.0;

        for(int st = 0; st < num_states; st++)
        {
            if( std::abs(kptr->Kstates[st].occupation[0]) < eps) 
            {
                vmax_st = st;
                break;
            }
        }

        for(int st = 0; st < num_states; st++)
        {
            if( std::abs(kptr->Kstates[st].occupation[0] - full_occ) > eps) 
            {
                cmin_st = st;
                break;
            }
        }

        //for(int stv = pct.gridpe; stv < vmax_st; stv+=pct.grid_npes)
        double tem = 0.0;
        for(int stv = pct.gridpe; stv < vmax_st; stv+=pct.grid_npes)
        {
            for(int stc = cmin_st; stc < num_states; stc++)
            {
                double ediff = kptr->Kstates[stc].eig[0] - kptr->Kstates[stv].eig[0];
                if(ediff > Emax || ediff < eps) continue;
                double occ_diff =  kptr->Kstates[stv].occupation[0] - kptr->Kstates[stc].occupation[0];
                int egrid_center = int(ediff/delta_e);
                for (int ie = - negrid; ie <= negrid; ie++)
                {
                    if (ie+egrid_center <0 || ie+egrid_center >= Epoints) continue;
                    double coeff = std::exp(-ie*ie * e2gaus) * occ_diff * kptr->kp.kweight;
                    for(int ix = 0; ix < 3; ix++)
                    {
                        for(int iy = 0; iy < 3; iy++)
                        {
                            tem = coeff * std::real(MyConj(Pmat[ix][stv * num_states + stc]) * Pmat[iy][stv * num_states+stc]);
                            epsilon[ix*3+iy][ie+egrid_center] += tem;
                        }
                    }
                }

                for (int ie = 0; ie < Epoints; ie++)
                {
                    double coeff = gaus_broad/(gaus_broad * gaus_broad+(ie-egrid_center) * (ie-egrid_center) * delta_e * delta_e) * occ_diff * kptr->kp.kweight/PI;
                    for(int ix = 0; ix < 3; ix++)
                    {
                        for(int iy = 0; iy < 3; iy++)
                        {
                            tem = coeff * std::real(MyConj(Pmat[ix][stv * num_states + stc]) * Pmat[iy][stv * num_states+stc]);
                            epsilon_lorentzian[ix*3+iy][ie] += tem;
                        }
                    }
                }

            }
        }


        //      for(int i = 0; i < 10; i++) 
        //      {
        //          rmg_printf("\n aaa" );
        //          for(int j = 0; j < 10; j++)
        //              rmg_printf(" %f ",std::real(Pxmat[i * num_states +j]));
        //      }
        //      for(int i = 0; i < 10; i++) 
        //      {
        //          rmg_printf("\n bbb");
        //          for(int j = 0; j < 10; j++)
        //              rmg_printf(" %f ",std::imag(Pxmat[i * num_states +j]));
        //      }
        if(pct.gridpe == 0) 
        {
            mkdir("Epsilon", S_IRWXU);
            int amode = S_IREAD | S_IWRITE;
            int kpt_glob = kpt + pct.kstart;
            size_t bytes;
            std::string filename = "Epsilon/EpsilonMat_spin"+std::to_string(pct.spinpe)+"_kpt"+std::to_string(kpt_glob);
            int fhand = open(filename.c_str(), O_CREAT | O_TRUNC | O_RDWR, amode);
            bytes = write(fhand, &num_states, sizeof(int));
            bytes = write(fhand, &kptr->kp.kweight, sizeof(double));
            for(int st = 0; st < num_states; st++) bytes = write(fhand, &kptr->Kstates[st].eig[0], sizeof(double));
            for(int st = 0; st < num_states; st++) bytes = write(fhand, &kptr->Kstates[st].occupation[0], sizeof(double));
            bytes = write(fhand, Pmat[0].data(), sizeof(KpointType) * num_states * num_states);
            bytes = write(fhand, Pmat[1].data(), sizeof(KpointType) * num_states * num_states);
            bytes = write(fhand, Pmat[2].data(), sizeof(KpointType) * num_states * num_states);
            if(bytes != sizeof(KpointType) * num_states * num_states)
            {
                rmg_error_handler (__FILE__, __LINE__, "size of writing epsilon mat is wrong \n");
            }
            close(fhand);
        }
    }
    for(int i = 0; i < 9; i++)
    {
        BlockAllreduce(epsilon[i].data(), (size_t)Epoints , pct.grid_comm);
        BlockAllreduce(epsilon[i].data(), (size_t)Epoints , pct.kpsub_comm);
        BlockAllreduce(epsilon_lorentzian[i].data(), (size_t)Epoints , pct.grid_comm);
        BlockAllreduce(epsilon_lorentzian[i].data(), (size_t)Epoints , pct.kpsub_comm);
    }

    double eps_tem[9];
    for(int ie =1; ie < Epoints; ie++)
    {

        for(int i = 0; i < 9; i++) eps_tem[i] = epsilon[i][ie];
        Rmg_Symm->symmetrize_tensor(eps_tem);

        double ene = ie * delta_e;
        for(int i = 0; i < 9; i++)
        {
            epsilon[i][ie] = eps_tem[i] * 4.0 * PI *PI/(ene * ene)/L->get_omega()/gaus_broad /std::sqrt(PI);

        }

        for(int i = 0; i < 9; i++) eps_tem[i] = epsilon_lorentzian[i][ie];
        Rmg_Symm->symmetrize_tensor(eps_tem);

        for(int i = 0; i < 9; i++)
        {
            epsilon_lorentzian[i][ie] = eps_tem[i] * 4.0 * PI *PI/(ene * ene)/L->get_omega();

        }
    }

    if(pct.gridpe == 0 && pct.kstart == 0) 
    {
        std::string filename = "Epsilon/epsilon_spin"+std::to_string(pct.spinpe)+".dat";
        FILE *eps_fi = fopen(filename.c_str(), "w");
        fprintf(eps_fi, "&& imag part epsilon2(omega) tensor, xx, yy,zz xy, xz, yz");
        for(int ie = 0; ie < Epoints; ie++)
        {
            fprintf(eps_fi, "\n%f    %e %e %e   %e %e %e", ie*delta_e*Ha_eV,epsilon[0][ie], epsilon[4][ie], epsilon[8][ie], epsilon[1][ie], epsilon[2][ie], epsilon[5][ie]);
        }
        
        // Kramers-Kronig to get the real part
        fprintf(eps_fi, "\n&& real part");
        for(int ie1 = 0; ie1 < Epoints; ie1++)
        {
            double ene1 = ie1 * delta_e;
            double tem = 0.0;

            for(int ie2 = 0; ie2 < Epoints; ie2++)
            {
                double ene2 = ie2 * delta_e;
                if(ie1 != ie2)
                    tem += (epsilon[0][ie2] + epsilon[1][ie2] + epsilon[3][ie2])/3.0 * ene2 * delta_e/(ene2* ene2 - ene1 * ene1);
            }
            double epsilon1 = 1.0 + 2.0/PI * tem;
            fprintf(eps_fi, "\n%f    %e ", ie1*delta_e*Ha_eV,epsilon1);

        }


        fclose(eps_fi);

        filename = "Epsilon/epsilonLorentzian_spin"+std::to_string(pct.spinpe)+".dat";
        eps_fi = fopen(filename.c_str(), "w");
        fprintf(eps_fi, "&& imag part epsilon2(omega) tensor, xx, yy,zz xy, xz, yz");
        for(int ie = 0; ie < Epoints; ie++)
        {
            fprintf(eps_fi, "\n%f    %e %e %e   %e %e %e", ie*delta_e*Ha_eV,epsilon_lorentzian[0][ie], epsilon_lorentzian[4][ie], epsilon_lorentzian[8][ie], epsilon_lorentzian[1][ie], epsilon_lorentzian[2][ie], epsilon_lorentzian[5][ie]);
        }
        fclose(eps_fi);
    }



    delete [] block_matrix;
}

