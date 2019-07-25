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


// 
//
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "main.h"
#include "common_prototypes.h"
int memory_steps;
double *Rold, *Fold, *direction, *alpha_lbfgs, *change_in_G, *change_in_R;
double *ro;
int reset_flag, fdstep, lineopt;
double finite_step,  maxmove, invcurv;

void lbfgs (double *posion, double *force, int num_ions, int num_images)
{


/**********************************************************************
!
! Limited-memory bfgs method, adapted from Vasp TST Tools 
!
!**********************************************************************
posion: cartesian coordinates of ions: 3 * num_ions *num_images;
force:  dimension of 3* num_ions * num_images, 
num_ions: total number of ions per image, if L-BFGS 
num_images: total number of images in NEB, if GL-BFGS 

*/

    double beta, a1,a2, curvature, fp1, fp2, Favg, step_size;
    int bound, im, jm;
    int item, i, ione =1;

    int tot_dim;

    double *R, *F;

    static int itr;

    tot_dim = 3 * num_ions * num_images;

    R=posion;
    F=force;


    if (fdstep)  // if no line minimizer is used (lineopt = 0), fdstep is alwalys 1
    {

        fdstep = 0;

        // check for reset of direction
        a1 = QMD_ddot (tot_dim, F, ione, Fold, ione);
        a1 = fabs(a1);

        a2 = QMD_ddot (tot_dim, Fold, ione, Fold, ione);
        a2 = fabs(a2);

        if (a2 < 1.0e-10) reset_flag = 1;
        if ((a1 > 0.5*a2) && lineopt) reset_flag = 1;


        if(reset_flag) 
        {

            itr=-1;
            reset_flag=0;
        }

        else

        {
            // find new direction
            itr = itr+1;
            if (itr < memory_steps) 

            {
                for( i = 0; i < tot_dim; i++)
                {
                    change_in_G[itr * tot_dim + i] = -(F[i] - Fold[i]); 
                    change_in_R[itr * tot_dim + i] = R[i] - Rold[i]; 
                }

                // apply periodic boundary conditions to these vector differences
                set_pbc(&change_in_R[itr * tot_dim],num_ions, num_images);



                a1 = QMD_ddot (tot_dim, &change_in_G[itr * tot_dim], ione, &change_in_R[itr * tot_dim], ione);
                ro[itr]=1.0/a1;
            }
            else
            {
                for(im = 0; im < memory_steps-1; im++)
                {
                    QMD_dcopy(tot_dim, &change_in_G[(im+1) * tot_dim], ione, &change_in_G[im * tot_dim], ione);
                    QMD_dcopy(tot_dim, &change_in_R[(im+1) * tot_dim], ione, &change_in_R[im * tot_dim], ione);
                    ro[im] = ro[im+1];
                }

                item = memory_steps -1;
                for( i = 0; i < tot_dim; i++)
                {
                    change_in_G[item * tot_dim + i] = -(F[i] - Fold[i]); 
                    change_in_R[item * tot_dim + i] = R[i] - Rold[i]; 
                }

                // apply periodic boundary conditions to these vector differences
                set_pbc(&change_in_R[item* tot_dim],num_ions, num_images);

                a1 = QMD_ddot (tot_dim, &change_in_G[item * tot_dim], ione, &change_in_R[item * tot_dim], ione);
                ro[item] = 1.0/a1;
            }
        }  //end if(reset_flag) 

        QMD_dcopy(tot_dim, R, ione, Rold, ione);
        QMD_dcopy(tot_dim, F, ione, Fold, ione);

        // compute Ho*g
        bound = memory_steps;
        if (itr < memory_steps) bound = itr+1;

        QMD_dcopy(tot_dim, F, ione, direction, ione);
        a1 = -1.0;
        QMD_dscal(tot_dim, a1, direction, ione);

        // calculate Ho*g
        for(im = 0; im < bound; im++)
        {
            jm = bound-im -1;
            alpha_lbfgs[jm] = ro[jm]*QMD_ddot(tot_dim, &change_in_R[jm* tot_dim], ione, direction, ione);
            a1 = - alpha_lbfgs[jm];
            QMD_daxpy(tot_dim, a1, &change_in_G[jm * tot_dim], ione, direction, ione);
        }


        a1 = invcurv;
        QMD_dscal(tot_dim, a1, direction, ione); //Ho=Ho*invcurv Ho is identity matrix


        for(im = 0; im < bound; im++)
        {
            beta = ro[im]*QMD_ddot(tot_dim, &change_in_G[im* tot_dim], ione, direction, ione);
            a1 =  alpha_lbfgs[im] - beta;
            QMD_daxpy(tot_dim, a1, &change_in_R[im * tot_dim], ione, direction, ione);
        }

        // direction down gradient
        a1 = -1.0;
        QMD_dscal(tot_dim, a1, direction, ione);

        if (lineopt) // line minimizer (using fdstep)
        {


            a1 = QMD_ddot (tot_dim, direction, ione, direction, ione);
            a1 = 1.0/sqrt(a1);

            QMD_dscal(tot_dim, a1, direction, ione); 


            // finite step down force
            a1 =  finite_step;
            QMD_daxpy(tot_dim, a1, direction, ione, R, ione);

        }
        else // use hessian to make step
        {
            a1 = QMD_ddot (tot_dim, direction, ione, direction, ione);
            step_size = fabs(sqrt(a1));
            if (step_size > sqrt(maxmove*maxmove*num_images)) 
            { 
                step_size=sqrt(maxmove*maxmove*num_images);
                a1 = step_size /sqrt(a1);
                QMD_dscal(tot_dim, a1, direction, ione); 
            }

            a1 = 1.0;
            QMD_daxpy(tot_dim, a1, direction, ione, R, ione);
            fdstep = 1;
        }
    }

    else  // Translation step for line-minimizer

    {
        fdstep = 1;

        // calculate curvature down direction
        fp1 = QMD_ddot (tot_dim, Fold, ione, direction, ione);
        fp2 = QMD_ddot (tot_dim, F, ione, direction, ione);
        curvature = (fp1-fp2)/finite_step;
        if (curvature < 0.0) 
        {
            step_size=maxmove;
        }
        else
        {
            Favg = 0.5*(fp1+fp2);
            step_size = Favg/curvature;
            if (fabs(step_size) > maxmove) 
            {
                step_size = (maxmove - finite_step) * fabs(step_size)/step_size;
            }
            else
            {
                step_size = step_size-0.5*finite_step ;
            }
        }


        // Move now from the configuration after the fd_step, so (*) has a "-" sign
        a1 = step_size;
        daxpy(&tot_dim, &a1, direction, &ione, R, &ione);

    }


}

