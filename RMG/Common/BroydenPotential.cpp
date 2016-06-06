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
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <boost/circular_buffer.hpp>



#define MAX_BROYDEN_ITER 12

void BroydenPotential(double *rho, double *new_rho, double *rhoc, double *vh_in, double *vh_out, int max_iter, bool reset)
{

   static boost::circular_buffer<double *> df(MAX_BROYDEN_ITER);
   static boost::circular_buffer<double *> dv(MAX_BROYDEN_ITER);
   static boost::circular_buffer<double *> dvh(MAX_BROYDEN_ITER);
   static int iter;
   double betamix[MAX_BROYDEN_ITER][MAX_BROYDEN_ITER];
   int pbasis = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
   int ld_betamix = MAX_BROYDEN_ITER;

   // Check if this is a reset request
   if(reset) {
       // Clear any old allocations
       while(df.size()) {
           delete [] df[0];
           df.pop_front();
       }
       while(dv.size()) {
           delete [] dv[0];
           dv.pop_front();
       }
       while(dvh.size()) {
           delete [] dvh[0];
           dvh.pop_front();
       }
       iter = 0;
       return;
   }


   // Set up arrays and get delta rho
   double *rhout = new double[pbasis];
   double *rhoin = new double[pbasis];
   for(int i = 0;i < pbasis;i++) rhoin[i] = rho[i];
   for(int i = 0;i < pbasis;i++) rhout[i] = new_rho[i];
   for(int i = 0;i < pbasis;i++) rhout[i] -= rhoin[i]; 


   // Check if it's time to remove old entries
   if(df.size() == (size_t)max_iter) {
       delete [] df[0];
       delete [] dv[0];
       delete [] dvh[0];
       df.pop_front();
       dv.pop_front();
       dvh.pop_front();
   }


   if(iter > 0) {
       for(int i = 0;i < pbasis;i++) df[df.size() - 1][i] -= rhout[i];
       for(int i = 0;i < pbasis;i++) dv[df.size() - 1][i] -= rhoin[i];
       for(int i = 0;i < pbasis;i++) dvh[df.size() - 1][i] -= (vh_out[i] - vh_in[i]);
   }

   int iter_used = df.size();

   // Create new entries
   double *df1 = new double[pbasis];
   for(int i = 0;i < pbasis;i++) df1[i] = rhout[i];
   df.push_back(df1);

   double *dv1 = new double[pbasis];
   for(int i = 0;i < pbasis;i++) dv1[i] = rhoin[i];
   dv.push_back(dv1);

   double *dvh1 = new double[pbasis];
   for(int i = 0;i < pbasis;i++) dvh1[i] = vh_out[i] - vh_in[i];
   dvh.push_back(dvh1);


   if(iter_used > 0) {

       for(int i = 0;i < iter_used;i++) {
           for(int j = 0;j < iter_used;j++) {
               betamix[j][i] = 0.0;
               for(int k = 0;k < pbasis;k++) betamix[j][i] += df[i][k] * dvh[j][k];
           }
       }

       MPI_Allreduce(MPI_IN_PLACE, betamix, MAX_BROYDEN_ITER*MAX_BROYDEN_ITER, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

       double work[MAX_BROYDEN_ITER*MAX_BROYDEN_ITER];
       int ipiv[MAX_BROYDEN_ITER*MAX_BROYDEN_ITER];
       int lwork = MAX_BROYDEN_ITER*MAX_BROYDEN_ITER;
       int info = 0;
       dgetrf(&iter_used, &iter_used, &betamix[0][0], &ld_betamix, ipiv,  &info );
       if(info)
           throw RmgFatalException() << "dgetrf failed " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

       dgetri(&iter_used, &betamix[0][0], &ld_betamix, ipiv, work, &lwork, &info );
       if(info)
           throw RmgFatalException() << "dgetri failed " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

       
       for(int i = 0;i < iter_used;i++) {
           work[i] = 0.0;
           for(int k = 0;k < pbasis;k++) work[i] += dvh[i][k] * rhout[k];
           MPI_Allreduce(MPI_IN_PLACE, &work[i], 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
       }

       for(int i = 0;i < iter_used;i++) {

           double gamma = 0.0;
           for(int j=0;j < iter_used;j++) gamma += work[j] * betamix[i][j];

           for(int k=0;k < pbasis;k++) rhout[k] -= gamma * df[i][k];
           for(int k=0;k < pbasis;k++) rhoin[k] -= gamma * dv[i][k];

       }

   }

   // New density. This is the place to do things with screening and frequency based mixing
   for(int k=0;k < pbasis;k++)  rhoin[k] =  rhoin[k] + ct.charge_broyden_scale*rhout[k];
   for(int k=0;k < pbasis;k++)  rho[k] = rhoin[k];

   delete [] rhoin;
   delete [] rhout;

   iter++;

}


