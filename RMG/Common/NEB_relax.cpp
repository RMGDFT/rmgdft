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


// Instantiate gamma and non-gamma versions
template void NEB_relax<double>(int , double *, double *, double *,
              double *, double *, double *, double *, Kpoint<double> **Kptr);
template void NEB_relax<std::complex<double> >(int , double *, double *, double *,
              double *, double *, double *, double *, Kpoint<std::complex<double> >**Kptr);



template <typename OrbitalType> void NEB_relax (int max_steps, double * vxc, double * vh, double * vnuc,
              double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr)
{

    int constrain;

    // atomic coordinates for left, self, and right images. 
    double *L_coor = new double[3*ct.num_ions];
    double *S_coor = new double[3*ct.num_ions];
    double *R_coor = new double[3*ct.num_ions];
    double *totale = new double[3*pct.images];
    double *all_frc = &totale[pct.images];
    double *path_length = &totale[2*pct.images];

    // total energy for left, self, and right images. 
    double L_total, S_total, R_total;

    double tmp_mag, max_frc, *fp;
    bool CONV_FORCE;

    MPI_Status status;
    if(pct.worldrank == 0) printf("\tEntering NEB routine.\n");

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

        if(pct.thisimg == 0 || pct.thisimg == pct.images -1)
            ct.constrainforces = 0;

        /* capture force constraint parameters from right and left data*/
        for(int img = 0; img < pct.images; img++) path_length[img] = 0.0;
        for (int count = 0; count < ct.num_ions; count++ )
        {
            double x = S_coor[3*count +0] - R_coor[3*count + 0];
            double y = S_coor[3*count +1] - R_coor[3*count + 1];
            double z = S_coor[3*count +2] - R_coor[3*count + 2];
            path_length[pct.thisimg] += x*x + y*y + z*z;
            /* put force constraints into control structure */
            Atoms[count].constraint.setA_weight = L_total;
            Atoms[count].constraint.setA_coord[0] = L_coor[3*count + 0]; 
            Atoms[count].constraint.setA_coord[1] = L_coor[3*count + 1]; 
            Atoms[count].constraint.setA_coord[2] = L_coor[3*count + 2]; 

            Atoms[count].constraint.setB_weight = R_total;
            Atoms[count].constraint.setB_coord[0] = R_coor[3*count + 0]; 
            Atoms[count].constraint.setB_coord[1] = R_coor[3*count + 1]; 
            Atoms[count].constraint.setB_coord[2] = R_coor[3*count + 2]; 

            /* zero velocities for every nudge */
            Atoms[count].velocity[0] = 0.0;
            Atoms[count].velocity[1] = 0.0;
            Atoms[count].velocity[2] = 0.0;
        }

        path_length[pct.thisimg] = sqrt(path_length[pct.thisimg]);
        /* Call fastrelax for max_md_steps steps */
        MPI_Barrier( MPI_COMM_WORLD );
        if(pct.worldrank == 0) rmg_printf("\tNEB call fast relax.\n");
        for (int count = 0; count < ct.num_ions; count++)
        {
            Atoms[count].constraint.forcemask[0] =0.0;
            Atoms[count].constraint.forcemask[1] =0.0;
            Atoms[count].constraint.forcemask[2] =0.0;
        }



        Relax (2, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);
        MPI_Barrier( MPI_COMM_WORLD );

        //rmg_lbfgs();

        /* Check for NEB convergence */
        /* Are we force converged? */
        tmp_mag = 0.0;
        max_frc = 0.0;
        for (int count = 0; count < ct.num_ions; count++)
        {
            double fx = Atoms[count].force[ct.fpt[0]][0] + Atoms[count].constraint.forcemask[0];
            double fy = Atoms[count].force[ct.fpt[0]][1] + Atoms[count].constraint.forcemask[1];
            double fz = Atoms[count].force[ct.fpt[0]][2] + Atoms[count].constraint.forcemask[2];
            tmp_mag =  fx*fx + fy*fy + fz*fz;
            if ( tmp_mag > max_frc )
                max_frc = tmp_mag;
        }

        max_frc = std::sqrt(max_frc);
        MPI_Barrier( MPI_COMM_WORLD );

        for(int pe = 0; pe < pct.images; pe++) totale[pe] = 0.0;
        for(int pe = 0; pe < pct.images; pe++) all_frc[pe] = 0.0;
        totale[pct.thisimg] = ct.TOTAL;
        all_frc[pct.thisimg] = max_frc;

        int count = 3 * pct.images;
        MPI_Allreduce(MPI_IN_PLACE, totale, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        MPI_Barrier( MPI_COMM_WORLD );
        if(pct.worldrank == 0)
        {
            printf("\n image     total energy(eV)     max_force \n");
            for(int img = 0; img < pct.images; img++)
            {
                printf("\n %d        %15.6e         %15.6e    %f", img, totale[img]/2, all_frc[img], path_length[img]);
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

    } 

    delete [] L_coor;
    delete [] S_coor;
    delete [] R_coor;
    delete [] totale;
}                               /* end neb_relax */


