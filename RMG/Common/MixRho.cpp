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
#include "RmgSumAll.h"
#include "blas.h"
#include "GlobalSums.h"
#include "Functional.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <boost/circular_buffer.hpp>



void mix_johnson(double *xm, double *fm, int NDIM, int ittot);



void MixRho (double * new_rho, double * rho, double *rhocore, double *vh_in, double *vh_out, double *rhoc, std::unordered_map<std::string, InputKey *>& ControlMap, bool reset)
{
    RmgTimer RT0("Mix rho");
    double t1, nspin = (ct.spin_flag + 1.0);
    int pbasis = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    int pbasis_noncoll = ct.noncoll_factor * ct.noncoll_factor * pbasis;

    //ct.mix = AutoMix (new_rho, rho, rhocore, vh_in, vh_out, rhoc, ControlMap, reset);

    if(Verify ("freeze_occupied", true, ControlMap)) return;

    /*Linear Mixing*/
    if (Verify("charge_mixing_type","Linear", ControlMap) || ct.charge_pulay_order == 1 || ((ct.scf_steps < ct.davidson_premg) && (ct.md_steps == 0) && (ct.runflag != RESTART )) || (ct.xc_is_hybrid && Functional::is_exx_active()))
    {
        if(reset) return;
        RmgTimer RT1("Mix rho: Linear");
        /* Scale old charge density first*/
        std::vector<double> drho;
        drho.resize(pbasis_noncoll);
        for(int ix = 0;ix < pbasis_noncoll;ix++) drho[ix] = new_rho[ix] - rho[ix];
        if(ct.drho_precond) Precond_drho(drho.data());
        t1 = ct.mix;
        int ione = 1;
        daxpy(&pbasis_noncoll, &t1, drho.data(), &ione, rho, &ione);

        rmg_printf("Charge density mixing: Linear with a constant of %.2f \n", t1);
    }
    else if (Verify("charge_mixing_type","Pulay", ControlMap))
    {
        RmgTimer RT1("Mix rho: Pulay");
        if(!Pulay_rho)
        {
            Pulay_rho = new PulayMixing(pbasis_noncoll, ct.charge_pulay_order, ct.charge_pulay_refresh,
                    ct.mix, ct.mix, pct.grid_comm);
            Pulay_rho->SetGspace(ct.drho_precond, ct.charge_pulay_Gspace, ct.drho_q0);

        }

        if(reset) {
            Pulay_rho->Refresh();
            return;
        }

        double mone = -1.0;
        int ione = 1;
        daxpy(&pbasis_noncoll, &mone, rho, &ione, new_rho, &ione);

        // rho_new store thr rho resudyke,
        if(ct.charge_pulay_Gspace) {
            Pulay_rho->Mixing_rhoG(rho, new_rho);
        }
        else {

            Pulay_rho->Mixing(rho, new_rho);
        }

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

    ct.min_rho = min;

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


// Function to automatically adjust mixing parameter
double AutoMix (double * new_rho, double * rho, double *rhocore, double *vh_in, double *vh_out, double *rhoc, std::unordered_map<std::string, InputKey *>& ControlMap, bool &reset)
{
    static std::vector<double> RMSdrho;
    static double adj_factor = 0.3;

    double drho = 0.0;
    int pbasis = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int pbasis_noncoll = ct.noncoll_factor * ct.noncoll_factor * pbasis;
    for(int ix=0;ix < pbasis_noncoll;ix++) drho += (new_rho[ix]-rho[ix])*(new_rho[ix]-rho[ix]);
    drho = RmgSumAll(drho, pct.grid_comm);
    drho = RmgSumAll(drho, pct.spin_comm);
    drho /= ct.psi_fnbasis;
    drho = sqrt(drho);
    RMSdrho.push_back(drho);
    if(ct.verbose && pct.gridpe == 0) printf("\nRMSdrho = %16.8e\n", drho);

    int j = RMSdrho.size();
    double newmix = ct.mix;
    if(j >= 5)
    {
        double ravg0 = (RMSdrho[j-4] + RMSdrho[j-3] + RMSdrho[j-2] + RMSdrho[j-1]) / 4.0;
        double ravg1 = (RMSdrho[j-5] + RMSdrho[j-4] + RMSdrho[j-3] + RMSdrho[j-2]) / 4.0;
        if(ravg0 > ravg1)
        {
            newmix = ct.mix*(1.0-adj_factor);
            RMSdrho.clear();
        }
        else
        {
//            newmix = std::min(1.0, ct.mix*((1.0-adj_factor)+adj_factor/2.0));
        }
        if(pct.gridpe == 0)
        {
            printf("\nOldmix = %7.4f, Newmix = %7.4f", ct.mix, newmix);
        }
    }
    return newmix;
}
