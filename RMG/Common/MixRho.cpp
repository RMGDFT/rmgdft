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

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cfloat>
#include <climits>
#include <unordered_map>
#include <typeinfo>
#include "const.h"
#include "InputKey.h"
#include "common_prototypes.h"
#include "RmgParallelFft.h"
#include "RmgException.h"
#include "transition.h"
#include "blas.h"
#include "GlobalSums.h"
#include "rmg_alloc.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <boost/circular_buffer.hpp>


static int pulay_step = 0;

void mix_johnson(double *xm, double *fm, int NDIM, int ittot);


void MixRho (double * new_rho, double * rho, double *rhocore, double *vh_in, double *vh_out, double *rhoc, std::unordered_map<std::string, InputKey *>& ControlMap, bool reset)
{
    RmgTimer RT0("Mix rho");
    double t1, nspin = (ct.spin_flag + 1.0);
    int pbasis = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    int length_x = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    int length_y = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    int length_z = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);


    if(Verify ("freeze_occupied", true, ControlMap)) return;

    /*Linear Mixing*/
    //if (Verify("charge_mixing_type","Linear", ControlMap) || ct.charge_pulay_order == 1 || (ct.rms > 5.0e-3) || (ct.scf_steps < 4))
    if (Verify("charge_mixing_type","Linear", ControlMap) || ct.charge_pulay_order == 1 || (ct.scf_steps < 4))
    {
        if(reset) return;
	RmgTimer RT1("Mix rho: Linear");
	/* Scale old charge density first*/
	t1 = ct.mix;
        for(int ix = 0;ix < pbasis;ix++) rho[ix] *= (1.0 - t1);

	/*Add the new density*/
        for(int ix = 0;ix < pbasis;ix++) rho[ix] += t1 * new_rho[ix];
        rmg_printf("Charge density mixing: Linear with a constant of %.2f \n", t1);
    }
    else if (Verify("charge_mixing_type","Pulay", ControlMap))
    {

        static double **rhohist=NULL, **residhist=NULL;
	RmgTimer RT1("Mix rho: Pulay");
        if (ct.charge_pulay_refresh)
            pulay_step = pulay_step % ct.charge_pulay_refresh;
        if(reset) {
            pulay_step = 0;
            return;
        }

        /*Use pulay mixing, result will be in rho*/
        pulay_rho(pulay_step, pbasis, length_x, length_y, length_z, 
                  new_rho, rho, ct.charge_pulay_order, &rhohist, &residhist, 
                  ct.charge_pulay_special_metrics, ct.charge_pulay_special_metrics_weight);
	    
            pulay_step++;
        rmg_printf("Charge density mixing: Pulay\n");

    }
    else if (Verify("charge_mixing_type","Broyden", ControlMap))
    {
	RmgTimer RT1("Mix rho: Broyden");
        BroydenPotential(rho, new_rho, rhoc, vh_in, vh_out, ct.charge_broyden_order, reset);
        if(reset) return;
        rmg_printf("Charge density mixing: Broyden\n");
    }
    

    /*Find charge minimum */
    double min = ZERO;
    double min2 = ZERO;
    for (int idx = 0; idx < pbasis; idx++)
    {
        if (rho[idx] < min)
            min = rho[idx];
        
	/*Here we check charge density with rhocore added*/
	if ((rho[idx] + rhocore[idx] / nspin) < min2)
         	min2 = rho[idx] + rhocore[idx] / nspin;
    }


    /*Find absolute minimum from all PEs */
    MPI_Allreduce (MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, pct.img_comm);
    MPI_Allreduce (MPI_IN_PLACE, &min2, 1, MPI_DOUBLE, MPI_MIN, pct.img_comm);

    if (min < ZERO)
    {
        rmg_printf ("\n charge density is NEGATIVE after interpolation, minimum is %e\n", min);
        rmg_printf (" minimum charge density with core charge added is %e\n", min2);
    }

}



/*
 *  The charge density mixing subroutine 
 *  The program is based on the mixing method by 
 *  D.D.Johnson, PRB38, 12807(1988). 
 *  
 *  Written by Wenchang Lu, Feb. 8, 2000 NCSU
 */

/*

                  mix_johnson.c

*/

/*  w1, w0,wj, mixing parameters                 */
#define         w1            0.5
#define         w0            0.3
#define         wj            1.0
#define       maxiter         7

double *un=NULL;
double *df=NULL;
double *betakn=NULL;
double *dff=NULL;
double *xm1=NULL;
double *fm1=NULL;
double *aij=NULL;
double *aij2=NULL;
 
void mix_johnson(double *xm, double *fm, int NDIM, int ittot)
{


    if(!un) un = new double[NDIM * maxiter]();
    if(!df) df = new double[NDIM * maxiter]();
    if(!betakn) betakn = new double[maxiter * maxiter]();
    if(!dff) dff = new double[NDIM * maxiter]();
    if(!xm1) xm1 = new double[NDIM];
    if(!fm1) fm1 = new double[NDIM];
    if(!aij) aij = new double[maxiter*maxiter];
    if(!aij2) aij2 = new double[maxiter*maxiter];

    double fnorm, aaa, tem;
    int  iter, i, j, idx;


    iter = ittot % maxiter;

    if(iter == 0)
    {

        for(i = 0; i < NDIM; i++) {
            xm1[i] = xm[i];
            fm1[i] = fm[i]; 
            xm [i] += w1*fm[i];
        }   /* endfor */

    }
    else
    {
        fnorm = 0;
        for(i = 0; i < NDIM; i++) {
            fm1[i] = fm[i] - fm1[i];
            fnorm += fm1[i]*fm1[i];
        } /* endfor */
        idx = 1;

        MPI_Allreduce(&fnorm, &tem, idx, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
        fnorm = tem;

        fnorm = sqrt(fnorm);

        for(i = 0; i < NDIM; i++) {
            dff[i + (iter-1) * NDIM] = fm1[i]/fnorm;
            un[i + (iter-1) * NDIM] = w1*dff[i + (iter-1) * NDIM] 
                + (xm[i] - xm1[i])/fnorm;

            xm1[i] = xm[i];
            fm1[i] = fm[i]; 
        }   /* endfor */

        /*
         * Calculate betakn   
         */

        if(iter == 1)
        {
            betakn[0] = 1.0/(w0*w0 + wj*wj);
        }
        else
        {

            for(i = 0; i < iter-1; i++){

                tem = 0;

                for(j = 0; j < NDIM; j++) {
                    tem += dff[j+i*NDIM]*dff[j+(iter-1)*NDIM];
                }   /* endfor  */    
                tem *= wj*wj;
                idx = 1;
                MPI_Allreduce(&tem, &aij[i], idx, MPI_DOUBLE,
                        MPI_SUM, MPI_COMM_WORLD);

            }   /* endfor */


            aaa =0.0;
            for(i = 0; i < iter-1; i++){

                tem = 0.0;
                for(j = 0; j < iter-1; j++) 
                    tem += betakn[i+ j*maxiter]*aij[j];
                aaa += aij[i]*tem;

            }  /* endfor  */

            betakn[iter-1 + (iter-1) *maxiter] = 1.0/(wj*wj-aaa);
            for(i = 0; i < iter-1; i++){

                aaa =0.0;
                for(j = 0; j < iter-1; j++) 
                    aaa += betakn[i+ j*maxiter]*aij[j]; 

                betakn[i + (iter-1) *maxiter] = -betakn[iter-1 + (iter-1) *maxiter]*aaa;
                betakn[iter-1 + i* maxiter] = betakn[i+(iter-1) *maxiter];

            }  /* endfor  */

            for(i = 0; i < iter-1; i++){
                for(j = 0; j < iter-1; j++)
                    betakn[i+ j* maxiter] += betakn[i+ (iter-1) *maxiter]
                        *betakn[iter-1+j*maxiter]/betakn[iter-1+ (iter-1) *maxiter];
            }  /* endfor */

        }

        /*
         * updating of betakn is finished
         */


        /*
         * Calculate the new charge density or vector "xm"
         */


        for(i = 0; i < iter; i++){
            tem = 0.0;
            for(j = 0; j < NDIM; j++) tem += dff[j+i*NDIM]*fm[j];
            idx = 1;
            MPI_Allreduce(&tem, &aij2[i], idx, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
        }  /* endfor */

        for(i = 0; i < iter; i++){
            tem = 0.0;
            for(j = 0; j <iter; j++) tem += betakn[i+j*maxiter]*aij2[j];
            aij[i] = tem;
        }  /* endfor */


        for(i = 0; i < NDIM; i++){
            aaa = 0.0;
            for(j = 0; j <iter; j++) aaa += un[i+j*NDIM]*aij[j];
            xm[i] += w1*fm[i] - aaa;
        }  /* endfor */

    }

}



