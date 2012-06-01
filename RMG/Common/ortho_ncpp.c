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

                    ortho_ncpp.c

   Performs orthogonalization and normalization of states for norm
   conserving pseudopotentials.


*/

#include "main.h"
#include <float.h>
#include <math.h>



void ortho_ncpp(STATE *states)
{

    int kpt;
    REAL time1;
    time1 = my_crtc();

    for(kpt = 0;kpt < ct.num_kpts;kpt++) {

        /* Do the orthogonalization here */  
        gram(&ct.kp[kpt], ct.vel, ct.num_states, ct.num_states, P0_BASIS, P0_BASIS);

    } /* end for */

    rmg_timings (ORTHO_TIME, (my_crtc () - time1));

} /* end ortho_full */

