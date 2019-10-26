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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "transition.h"
#include "Atomic.h"
#include "RmgParallelFft.h"
#include "prototypes_rmg.h"
#include "GlobalSums.h"
#include "Neb.h"

template Neb<double>::Neb(BaseGrid &BG, int num_images, int max_steps, std::string input_initial, 
        std::string input_final, double totale_initial, double totale_final);
template Neb<std::complex<double>>::Neb(BaseGrid &BG, int num_images, int max_steps, std::string input_initial, 
        std::string input_final, double totale_initial, double totale_final);
template <class T> Neb<T>::Neb( BaseGrid &G_in, int num_images, int max_steps, std::string input_initial, 
        std::string input_final, double totale_initial, double totale_final):BG(G_in)
{

    // num_images: exclude inital and final images.
    this->num_images = num_images;
    this->max_steps = max_steps;
    this->totale_initial = totale_initial;
    this->totale_final = totale_final;

    L_coor = new double[3*ct.num_ions];
    S_coor = new double[3*ct.num_ions];
    R_coor = new double[3*ct.num_ions];
    totale = new double[this->num_images+2];
    all_frc = new double[this->num_images+2];
    path_length = new double[this->num_images +2];


    std::set<std::string> SpeciesTypes;
    std::list<std::string> IonSpecies;

    // Read inital image atoms
    ReadRmgAtoms(const_cast<char*>(input_initial.c_str()), SpeciesTypes, IonSpecies, Atoms_initial, ct);
    ReadRmgAtoms(const_cast<char*>(input_final.c_str()), SpeciesTypes, IonSpecies, Atoms_final, ct);

}

template Neb<double>::~Neb();
template Neb<std::complex<double>>::~Neb();
template <class T> Neb<T>::~Neb()
{
    delete [] L_coor;
    delete [] S_coor;
    delete [] R_coor;
    delete [] totale;
    delete [] all_frc;
    delete [] path_length;

}

template void Neb<double>::relax(double *, double *, double *,
        double *, double *, double *, double *, Kpoint<double> **Kptr);
template void Neb<std::complex<double>>::relax(double *, double *, double *,
        double *, double *, double *, double *, Kpoint<std::complex<double>> **Kptr);

template <class T> void Neb<T>::relax (double * vxc, double * vh, double * vnuc,
        double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<T> **Kptr)
{

    // total energy for left, self, and right images. 
    double L_total, S_total, R_total;

    double tmp_mag, max_frc, *fp;
    bool CONV_FORCE;

    MPI_Status status;
    if(pct.worldrank == 0) printf("\tEntering NEB routine.\n");

    ct.constrainforces = 5;
    /* Loop NEB relaxations */
    for (int neb_step = 0; neb_step < max_steps; neb_step++)
    {
        if(pct.worldrank == 0) printf ("\nNEBrlx step: ----------  %d  ----------\n", neb_step);

        /* pack coordinates for mpi transfer */
        for (int count = 0; count < ct.num_ions; count++ )
        {
            S_coor[3*count + 0] = Atoms[count].crds[0];
            S_coor[3*count + 1] = Atoms[count].crds[1];
            S_coor[3*count + 2] = Atoms[count].crds[2];
        }

        // communicate with left and right images
        int num_coor = 3 * Atoms.size();
        int tag = pct.gridpe;
        MPI_Sendrecv(S_coor, num_coor, MPI_DOUBLE, pct.left_img_rank, tag,
                R_coor, num_coor, MPI_DOUBLE, pct.right_img_rank, tag, MPI_COMM_WORLD, &status);

        MPI_Sendrecv(S_coor, num_coor, MPI_DOUBLE, pct.right_img_rank, tag,
                L_coor, num_coor, MPI_DOUBLE, pct.left_img_rank, tag, MPI_COMM_WORLD, &status);

        S_total = ct.TOTAL;

        MPI_Sendrecv(&S_total, 1, MPI_DOUBLE, pct.left_img_rank, tag,
                &R_total, 1, MPI_DOUBLE, pct.right_img_rank, tag, MPI_COMM_WORLD, &status);

        MPI_Sendrecv(&S_total, 1, MPI_DOUBLE, pct.right_img_rank, tag,
                &L_total, 1, MPI_DOUBLE, pct.left_img_rank, tag, MPI_COMM_WORLD, &status);

        // for the first image, its left is the initial image
        if(pct.thisimg == 0)   
        {
            L_total = this->totale_initial;
            for (int count = 0; count < ct.num_ions; count++ )
            {
                L_coor[3*count + 0] = Atoms_initial[count].crds[0];
                L_coor[3*count + 1] = Atoms_initial[count].crds[1];
                L_coor[3*count + 2] = Atoms_initial[count].crds[2];
            }

        }

        // for the last image, its right is the final image
        if(pct.thisimg == num_images-1) 
        {
            R_total = this->totale_final;
            for (int count = 0; count < ct.num_ions; count++ )
            {
                R_coor[3*count + 0] = Atoms_final[count].crds[0];
                R_coor[3*count + 1] = Atoms_final[count].crds[1];
                R_coor[3*count + 2] = Atoms_final[count].crds[2];
            }
        }

        /* capture force constraint parameters from right and left data*/
        for(int img = 0; img < pct.images +2; img++) path_length[img] = 0.0;

        double path_length0 = 0;
        for (int count = 0; count < ct.num_ions; count++ )
        {
            double rdiff[3], rdiff_crys[3];
            rdiff[0] = R_coor[3*count +0] - S_coor[3*count + 0];
            rdiff[1] = R_coor[3*count +1] - S_coor[3*count + 1];
            rdiff[2] = R_coor[3*count +2] - S_coor[3*count + 2];

            Rmg_L.to_crystal_half(rdiff_crys, rdiff);
            Rmg_L.to_cartesian(rdiff_crys, rdiff);

            path_length[pct.thisimg+1] += rdiff[0] * rdiff[0] + rdiff[1] * rdiff[1] + rdiff[2] * rdiff[2]; 

            Atoms[count].constraint.setB_weight = R_total;
            Atoms[count].constraint.setB_coord[0] = S_coor[3*count + 0] + rdiff[0]; 
            Atoms[count].constraint.setB_coord[1] = S_coor[3*count + 1] + rdiff[1]; 
            Atoms[count].constraint.setB_coord[2] = S_coor[3*count + 2] + rdiff[2]; 
            /* put force constraints into control structure */
            rdiff[0] = L_coor[3*count +0] - S_coor[3*count + 0];
            rdiff[1] = L_coor[3*count +1] - S_coor[3*count + 1];
            rdiff[2] = L_coor[3*count +2] - S_coor[3*count + 2];

            Rmg_L.to_crystal_half(rdiff_crys, rdiff);
            Rmg_L.to_cartesian(rdiff_crys, rdiff);
            if(pct.thisimg == 0) path_length[0] += 
                rdiff[0] * rdiff[0] + rdiff[1] * rdiff[1] + rdiff[2] * rdiff[2]; 
            Atoms[count].constraint.setA_weight = L_total;
            Atoms[count].constraint.setA_coord[0] = S_coor[3*count + 0] + rdiff[0]; 
            Atoms[count].constraint.setA_coord[1] = S_coor[3*count + 1] + rdiff[1]; 
            Atoms[count].constraint.setA_coord[2] = S_coor[3*count + 2] + rdiff[2]; 


        }

        path_length[pct.thisimg+1] = sqrt(path_length[pct.thisimg+1]);
        path_length[0] = sqrt(path_length[0]);
        /* Call fastrelax for max_md_steps steps */
        MPI_Barrier( MPI_COMM_WORLD );
        if(pct.worldrank == 0) rmg_printf("\tNEB call fast relax.\n");
        for (int count = 0; count < ct.num_ions; count++)
        {
            Atoms[count].constraint.forcemask[0] =0.0;
            Atoms[count].constraint.forcemask[1] =0.0;
            Atoms[count].constraint.forcemask[2] =0.0;
        }


        static double *rhodiff;

        Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr, true);
        WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);

        tmp_mag = 0.0;
        max_frc = 0.0;
        for (int count = 0; count < ct.num_ions; count++)
        {
            double fx = Atoms[count].force[ct.fpt[0]][0];
            double fy = Atoms[count].force[ct.fpt[0]][1];
            double fz = Atoms[count].force[ct.fpt[0]][2];
            tmp_mag =  fx*fx + fy*fy + fz*fz;
            if ( tmp_mag > max_frc )
                max_frc = tmp_mag;
        }

        max_frc = std::sqrt(max_frc);
        MPI_Barrier( MPI_COMM_WORLD );

        for(int pe = 0; pe < pct.images+2; pe++) totale[pe] = 0.0;
        for(int pe = 0; pe < pct.images+2; pe++) all_frc[pe] = 0.0;
        totale[pct.thisimg + 1] = ct.TOTAL;
        all_frc[pct.thisimg + 1] = max_frc;

        int count = pct.images + 2;
        if(pct.img_cross_comm == MPI_COMM_NULL) 
        {
            printf("\n NEB runs without setting img_cross_comm \n");
            printf("\n either num_processors for images are not equal or ct.image_per_node !=1\n");
            exit(0);
        }
        MPI_Allreduce(MPI_IN_PLACE, totale, count, MPI_DOUBLE, MPI_SUM, pct.img_cross_comm);
        MPI_Allreduce(MPI_IN_PLACE, all_frc, count, MPI_DOUBLE, MPI_SUM, pct.img_cross_comm);
        MPI_Allreduce(MPI_IN_PLACE, path_length, count, MPI_DOUBLE, MPI_SUM, pct.img_cross_comm);
        totale[0] = this->totale_initial;
        totale[this->num_images+1] = this->totale_final;

        MPI_Barrier( MPI_COMM_WORLD );
        if(pct.worldrank == 0)
        {
            double max_tote = *std::max_element(totale, totale+(pct.images+2));
            printf("\n Total Energy of Initial:  %f Ha", totale[0]);
            printf("\n Barrier from  left:  %f eV", (max_tote - totale[0])*Ha_eV);
            printf("\n Barrier from right:  %f eV", (max_tote - totale[(pct.images+2)-1])*Ha_eV);

            printf("\n image     total energy(eV)     max_force  path_length (au)\n");
            for(int img = 0; img < pct.images+2; img++)
            {
                printf("\n %d        %15.6e    %15.6e    %f", img, (totale[img]-totale[0])*Ha_eV, all_frc[img], path_length[img]);
            }
        }

        for(int img = 0; img < pct.images; img++)
        {
            if(max_frc < all_frc[img]) max_frc = all_frc[img];
        }

        CONV_FORCE = (max_frc < ct.thr_frc);

        if(CONV_FORCE) 
        {
            if(pct.worldrank == 0) printf("\n NEB converged\n");
            break;
        }



        // Get atomic rho for this ionic configuration and subtract from current rho
        int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
        double *arho = new double[FP0_BASIS];
        LcaoGetAtomicRho(arho);

        // If first step allocate rhodiff
        for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] -= arho[idx];

        if(rhodiff == NULL)
        {
            rhodiff = new double[FP0_BASIS];
            for(int idx = 0;idx < FP0_BASIS;idx++) rhodiff[idx] = rho[idx];
        }
        else
        {
            double *trho = new double[FP0_BASIS];
            for(int idx = 0;idx < FP0_BASIS;idx++) trho[idx] = rho[idx];
            for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] = 2.0*rho[idx] - rhodiff[idx];
            for(int idx = 0;idx < FP0_BASIS;idx++) rhodiff[idx] = trho[idx];
            delete [] trho;
        }

        /* not done yet ? => move atoms */
        /* move the ions */

        switch(ct.relax_method)
        {

            case FASTRELAX:
                fastrelax (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                break;
            case FIRE:
                fire (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                break;
            case QUICK_MIN:
                fastrelax (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                break;
            case MD_MIN:
                fastrelax (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                break;
            case LBFGS:
                rmg_lbfgs();
                break;
            default:
                rmg_error_handler (__FILE__, __LINE__, "Undefined MD method");
        }

        /* Update items that change when the ionic coordinates change */
        RmgTimer *RT0=new RmgTimer("1-TOTAL: run: ReinitIonicPotentials");
        ReinitIonicPotentials (Kptr, vnuc, rhocore, rhoc);
        delete RT0;

        // Reset mixing
        MixRho(NULL, NULL, NULL, NULL, NULL, NULL, Kptr[0]->ControlMap, true);

        // Get atomic rho for new configuration and add back to rho
        LcaoGetAtomicRho(arho);
        for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] += arho[idx];
        delete [] arho;

        // Extrapolate orbitals after first step
        ExtrapolateOrbitals(ct.outfile, Kptr);

    }

}                               /* end neb_relax */

