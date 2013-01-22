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
extern REAL *global_matrix;


void gram(KPOINT *kpt, REAL vel, int numst, int maxst, int numpt, int maxpt)
{

   int st, st1, info, length, idx, idj, omp_tid;
   int ione = 1, block, num_blocks, block_pts;
   REAL alpha = -1.0;
   REAL zero = 0.0;
   REAL one = 1.0;
   REAL tmp, *c, *tarr, *darr, *sarr;
   REAL time1;
   STATE *sp;
   char *transn = "n";
   char *transt = "t";
   char *uplo = "l";
   char *uphi = "u";
   int pbasis =pct.P0_BASIS;


   sp = kpt->kstate;
   c = sp->psiR;

   my_malloc(tarr, numst , REAL);

#if MD_TIMERS
   time1 = my_crtc();
#endif
   ssyrk( uplo, transt, &numst, &numpt, &one, c, &numpt, 
               &zero, global_matrix, &numst);
#if MD_TIMERS
   rmg_timings (ORTHO_GET_OVERLAPS, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   /* get the global part */
   length = maxst * maxst;
   global_sums(global_matrix, &length, pct.grid_comm);

#if MD_TIMERS
   rmg_timings (ORTHO_GLOB_SUM, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   /* compute the cholesky factor of the overlap matrix */
   cholesky(global_matrix, maxst);

#if MD_TIMERS
   rmg_timings (ORTHO_CHOLESKY, (my_crtc () - time1));
   time1 = my_crtc();
#endif


   // Get inverse of diagonal elements
   for(st = 0;st < numst;st++) tarr[st] = 1.0 / global_matrix[st + maxst * st];

// This code may look crazy but there is a method to the madness. We copy a slice
// of the wavefunction array consisting of the values for all orbitals of a given
// basis point into a temporary array. Then we do the updates on each slice and
// parallelize over slices with OpenMP. This produces good cache behavior
// and excellent parformance on the XK6.

#pragma omp parallel private(idx,st,st1,omp_tid,sarr)
{
       omp_tid = omp_get_thread_num();
       if(omp_tid == 0) my_malloc(darr, numst * omp_get_num_threads(), REAL);
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

#if MD_TIMERS
   rmg_timings (ORTHO_UPDATE_WAVES, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   tmp = 1.0 / sqrt(vel);
   idx = numst * maxpt;
   QMD_sscal(idx, tmp, c, ione);

#if MD_TIMERS
   rmg_timings (ORTHO_NORM_PSI, (my_crtc () - time1));
   time1 = my_crtc();
#endif


   my_free(darr);
   my_free(tarr);

} /* end of gram */


