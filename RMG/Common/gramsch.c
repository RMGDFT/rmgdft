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





void gram(KPOINT *kpt, REAL vel, int numst, int maxst, int numpt, int maxpt)
{

   int st, st1, info, length, idx, idj;
   int ione = 1;
   REAL alpha = -1.0;
   REAL zero = 0.0;
   REAL one = 1.0;
   REAL tmp, *darr, *c, *tarr;
   REAL time1;
   STATE *sp;
   char *transn = "n";
   char *transt = "t";
   char *uplo = "l";
   char *uphi = "u";
   int pbasis = P0_BASIS;


   sp = kpt->kstate;
   c = sp->psiR;

   my_malloc(darr, ct.num_states * ct.num_states, REAL);

#if MD_TIMERS
   time1 = my_crtc();
#endif
   ssyrk( uplo, transt, &numst, &numpt, &one, c, &numpt, 
               &zero, darr, &numst);
#if MD_TIMERS
   rmg_timings (ORTHO_GET_OVERLAPS, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   /* get the global part */
   length = maxst * maxst;
   global_sums(darr, &length, pct.grid_comm);

#if MD_TIMERS
   rmg_timings (ORTHO_GLOB_SUM, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   /* compute the cholesky factor of the overlap matrix */
   cholesky(darr, maxst);

#if MD_TIMERS
   rmg_timings (ORTHO_CHOLESKY, (my_crtc () - time1));
   time1 = my_crtc();
#endif

   /* apply inverse of cholesky factor to states */
   for (st = 0; st < numst; st++) {

      /* normalize c[st] */
      tmp = 1.0 / darr[st + maxst * st];
      QMD_sscal( numpt, tmp, &c[st * maxpt], ione);

      /* subtract the projection along c[st] from the remaining vectors */
      idx = numst - st - 1;
      if(idx)
          sger(&numpt, &idx, &alpha, &c[st * maxpt], &ione, 
               &darr[(st+1) + maxst*st], &ione, &c[(st+1) * maxpt], &maxpt);
 
   } /* end of for */

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

} /* end of gram */


