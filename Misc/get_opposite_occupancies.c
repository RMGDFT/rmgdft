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

/****f* QMD-MGDFT/scf.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   Obtaines occupations for wavefunctions with opposite spin
 * INPUTS
 *   states: points to orbital structure
 *   vxc: exchange correlation potential
 *   vh:  Hartree potential
 *   vnuc: pseudopotential
 *   rho: total valence charge density
 *   rhocore:  core charge density
 *   rhoc: Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 *   CONVERGENCE: 1 converged, 0 not
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c quench.c
 * CHILDREN
 *   
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void get_opposite_occupancies (STATE * states)
{

    MPI_Status status; 
    double *occ_sd, *occ_rv;   
    int st, st1, kpt;


    /* allocate memory for occupation send array and receive array */
    my_malloc (occ_sd, 2 * ct.num_kpts * ct.num_states, double);
    occ_rv = occ_sd + ct.num_kpts * ct.num_states;


    /*Prepare the sending buffer of occupancies */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            occ_sd[st] = ct.kp[pct.kstart + kpt].kstate[st1].occupation[0];
            st += 1;
        }	
    }


    /*Communicate for spin up and spin down occupations*/    
    MPI_Sendrecv(occ_sd, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
            occ_rv, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);


    /* Unpack the received occupations to state structure */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            ct.kp[pct.kstart + kpt].kstate[st1].occupation[1] = occ_rv[st];
            st += 1;

        }	
    } 


    my_free (occ_sd);
}                               /* end scf */


/******/
