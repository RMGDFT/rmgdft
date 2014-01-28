/*
 * QMD  Quantum molecular dynamics package.
 * Version: 2.1.3
 * Copyright (C) 1995  Emil Briggs
 * Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                     Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */


/*

                    gramsch.c


   Performs gram-schmidt orthogonalization of the wavefunctions



*/


#include "main.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>


// Also used by subdiag but never at the same time.
extern rmg_double_t *global_matrix;


void gram(KPOINT *kpt, rmg_double_t vel, int numst, int maxst, int numpt, int maxpt)
{

   int st, st1, info, length, idx, idj, omp_tid;
   int ione = 1, block, num_blocks, block_pts;
   rmg_double_t alpha = -1.0;
   rmg_double_t zero = 0.0;
   rmg_double_t one = 1.0;
   rmg_double_t tmp, *c, *tarr, *darr, *sarr;
   rmg_double_t time1;
   STATE *sp;
   char *transn = "n";
   char *transt = "t";
   char *uplo = "l";
   char *uphi = "u";
   int pbasis =pct.P0_BASIS;
#if GPU_ENABLED
   cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
   cublasFillMode_t cuplo=CUBLAS_FILL_MODE_LOWER;
#endif


   sp = kpt->kstate;
   c = sp->psiR;

   my_malloc(tarr, numst , rmg_double_t);

#if MD_TIMERS
   time1 = my_crtc();
#endif
#if GPU_ENABLED
   cublasSetVector( pbasis * numst, sizeof( rmg_double_t ), c, ione, ct.gpu_states, ione );
   cublasDsyrk_v2 (ct.cublas_handle, cuplo, cu_transT, numst, numpt, &one, ct.gpu_states, numpt, &zero, ct.gpu_global_matrix, numst);
   cublasGetVector( numst * numst, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );
#else
   ssyrk( uplo, transt, &numst, &numpt, &one, c, &numpt, 
               &zero, global_matrix, &numst);
#endif

#if MD_TIMERS
   rmg_timings (ORTHO_GET_OVERLAPS, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   /* get the global part */
   length = maxst * maxst;
   MPI_Allreduce(MPI_IN_PLACE, global_matrix, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);


#if MD_TIMERS
   rmg_timings (ORTHO_GLOB_SUM, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   /* compute the cholesky factor of the overlap matrix */
#if GPU_ENABLED && MAGMA_LIBS
   cublasSetVector( numst * numst, sizeof( rmg_double_t ), global_matrix, ione, ct.gpu_global_matrix, ione );
   magma_dpotrf_gpu('L', numst, ct.gpu_global_matrix, numst, &info);
   cublasGetVector( numst * numst, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );
#else
   cholesky(global_matrix, maxst);
#endif

#if MD_TIMERS
   rmg_timings (ORTHO_CHOLESKY, (my_crtc () - time1));
   time1 = my_crtc();
#endif


   // Get inverse of diagonal elements
   for(st = 0;st < numst;st++) tarr[st] = 1.0 / global_matrix[st + maxst * st];

#if GPU_ENABLED

   /* apply inverse of cholesky factor to states */
   for (st = 0; st < numst; st++) {

      /* normalize c[st] */
//      QMD_dscal( numpt, tarr[st], &c[st * maxpt], ione);
      cublasDscal_v2(ct.cublas_handle, numpt, &tarr[st], &ct.gpu_states[st * maxpt], ione);

      /* subtract the projection along c[st] from the remaining vectors */
      idx = numst - st - 1;
      if(idx) {
          cublasDger_v2 (ct.cublas_handle, numpt, idx, &alpha, &ct.gpu_states[st * maxpt], ione,
                         &ct.gpu_global_matrix[(st+1) + maxst*st], ione, &ct.gpu_states[(st+1) * maxpt], maxpt);
//          dger_(&numpt, &idx, &alpha, &c[st * maxpt], &ione,
//               &global_matrix[(st+1) + maxst*st], &ione, &c[(st+1) * maxpt], &maxpt);
      }

   } /* end of for */

   cublasGetVector( pbasis * numst, sizeof( rmg_double_t ), ct.gpu_states, ione, c, ione );
#else

// This code may look crazy but there is a method to the madness. We copy a slice
// of the wavefunction array consisting of the values for all orbitals of a given
// basis point into a temporary array. Then we do the updates on each slice and
// parallelize over slices with OpenMP. This produces good cache behavior
// and excellent parformance on the XK6.

#pragma omp parallel private(idx,st,st1,omp_tid,sarr)
{
       omp_tid = omp_get_thread_num();
       if(omp_tid == 0) my_malloc(darr, numst * omp_get_num_threads(), rmg_double_t);
#pragma omp barrier
       
#pragma omp for schedule(static, 1) nowait
   for(idx = 0;idx <pct.P0_BASIS;idx++) {

       sarr = &darr[omp_tid*numst];

       for (st = 0; st < numst; st++) sarr[st] = c[st*maxpt + idx];

       for (st = 0; st < numst; st++) {

           sarr[st] *= tarr[st];

           for (st1 = st+1; st1 < numst; st1++) {
               sarr[st1] -= global_matrix[st1 + maxst*st] * sarr[st];
           }

       }

       for (st = 0; st < numst; st++) c[st*maxpt + idx] = sarr[st];

   }
}
   my_free(darr);
#endif

#if MD_TIMERS
   rmg_timings (ORTHO_UPDATE_WAVES, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   tmp = 1.0 / sqrt(vel);
   idx = numst * maxpt;
   QMD_dscal(idx, tmp, c, ione);

#if MD_TIMERS
   rmg_timings (ORTHO_NORM_PSI, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   my_free(tarr);

} /* end of gram */


