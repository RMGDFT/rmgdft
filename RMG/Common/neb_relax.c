/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/neb_relax.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 2001  Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void neb_relax (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
 *                REAL *rho, REAL *rhocore, REAL *rhoc)
 *   drive routine for elementary nudged elastic band relax.
 * INPUTS
 *   states: all wave functions (see main.h)
 *   vxc: exchange-correlation potential
 *   vh:  Hartree potential
 *   vnuc: pseudopotential
 *   rho:  total charge density
 *   rhocore: charge density of core electrons, only useful when we 
 *            include non-linear core correction for pseudopotential.
 *   rhoc:    compensating charge density
 * OUTPUT
 *   all the above inputs are updated
 * PARENTS
 *   main.c
 * CHILDREN
 *   quench.c to_crystal.c to_cartesian.c init_nuc.c get_nlop.c scf.c
 *   get_te.c force.c rmg_fastrelax.c
 * SOURCE
 */

#include "main.h"
#include <math.h>
#include <stdbool.h>

#define LEFT 0
#define SELF 1
#define RIGHT 2

#define X 0
#define Y 1
#define Z 2

void neb_relax (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc,
              REAL * rho, REAL * rho_oppo, REAL * rhocore, REAL * rhoc)
{
    /* This may need to be malloced if we start using large ion counts */
    int constrain, count, neb_steps = 0, img_rank_map[3];
    REAL imgA[3*ct.num_ions], imgB[3*ct.num_ions], imgC[3*ct.num_ions];
    REAL tmp_mag, max_frc, *fp, *L_ptr, *S_ptr, *R_ptr;
	REAL L_total, S_total, R_total;
    bool CONV_FORCE, DONE = false;
    MPI_Request req[3];
    MPI_Status status;
	printf("\tEntering NEB routine.\n");

    if ( pct.imgpe == 0 )
    {
        MPI_Cart_shift (pct.img_topo_comm, 0, 1, &img_rank_map[LEFT], &img_rank_map[RIGHT] );
        MPI_Comm_rank (pct.img_topo_comm, &img_rank_map[SELF]);
    }
	printf("\tNEB routine: got adjacent image ranks in topo_comm.\n");

    /* Loop NEB relaxations */
    while ( !DONE )
    {
		printf ("\nNEBrlx: ---------- [md: %d/%d] ----------\n", ct.md_steps/ct.max_md_steps, ct.max_rmg_steps);
        L_ptr = imgA;
        S_ptr = imgB;
        R_ptr = imgC;

        /* get adjacent image ion positions */
        if ( pct.imgpe == 0 )
        {
			printf("\tNEB adjacent image coords L:%d, S:%d, R:%d.\n", img_rank_map[LEFT], img_rank_map[SELF], img_rank_map[RIGHT]);
			fflush (NULL); fsync( fileno(ct.logfile) );
            /* pack coordinates for mpi transfer */
            for ( count = 0; count < ct.num_ions; count++ )
            {
                S_ptr[3*count + X] = ct.ions[count].crds[X];
                S_ptr[3*count + Y] = ct.ions[count].crds[Y];
                S_ptr[3*count + Z] = ct.ions[count].crds[Z];
            }

            if ( img_rank_map[LEFT] != MPI_PROC_NULL  ) 
			{
				MPI_Isend ( S_ptr, 3*ct.num_ions, MPI_DOUBLE, img_rank_map[LEFT], LEFT, pct.img_topo_comm, &req[SELF] );
				MPI_Irecv ( L_ptr, 3*ct.num_ions, MPI_DOUBLE, img_rank_map[LEFT], RIGHT, pct.img_topo_comm, &req[LEFT] );
			}

            if ( img_rank_map[RIGHT] != MPI_PROC_NULL  ) 
			{
				MPI_Isend ( S_ptr, 3*ct.num_ions, MPI_DOUBLE, img_rank_map[RIGHT], RIGHT, pct.img_topo_comm, &req[SELF] );
				MPI_Irecv ( R_ptr, 3*ct.num_ions, MPI_DOUBLE, img_rank_map[RIGHT], LEFT, pct.img_topo_comm, &req[RIGHT] );
			}

			/* if data from left,wait; else use self data as left data */
            if ( img_rank_map[LEFT] != MPI_PROC_NULL  ) 
			{
				MPI_Wait( &req[LEFT], &status );
			}
			else
            {
                L_ptr = imgB;
                S_ptr = imgA;
            }

			/* if data from right,wait; else use self data as right data */
            if ( img_rank_map[RIGHT] != MPI_PROC_NULL  ) 
			{
				MPI_Wait( &req[RIGHT], &status );
			}
			else
            {
                R_ptr = imgB;
                S_ptr = imgC;
            }


			/* Set total energy for transfer */
			S_total = ct.TOTAL;

            if ( img_rank_map[LEFT] != MPI_PROC_NULL  ) 
			{
				MPI_Isend ( &S_total, 1, MPI_DOUBLE, img_rank_map[LEFT], LEFT, pct.img_topo_comm, &req[SELF] );
				MPI_Irecv ( &L_total, 1, MPI_DOUBLE, img_rank_map[LEFT], RIGHT, pct.img_topo_comm, &req[LEFT] );
			}

            if ( img_rank_map[RIGHT] != MPI_PROC_NULL  ) 
			{
				MPI_Isend ( &S_total, 1, MPI_DOUBLE, img_rank_map[RIGHT], RIGHT, pct.img_topo_comm, &req[SELF] );
				MPI_Irecv ( &R_total, 1, MPI_DOUBLE, img_rank_map[RIGHT], LEFT, pct.img_topo_comm, &req[RIGHT] );
			}

			/* if data from left,wait; else use self data as left data */
            if ( img_rank_map[LEFT] != MPI_PROC_NULL  ) 
			{
				MPI_Wait( &req[LEFT], &status );
			}
			else
            {
                L_total = S_total;
            }

			/* if data from right,wait; else use self data as right data */
            if ( img_rank_map[RIGHT] != MPI_PROC_NULL  ) 
			{
				MPI_Wait( &req[RIGHT], &status );
			}
			else
            {
                R_total = S_total;
            }

            /* End images are allowed to relax without constraint */
            if ( img_rank_map[RIGHT] == MPI_PROC_NULL || img_rank_map[LEFT] == MPI_PROC_NULL  ) 
			{
                constrain = 0;
            }
            else
            {
                constrain = ct.constrainforces;
            }

        }
		/* broadcast force constraint parameters to image procs */
		MPI_Bcast( L_ptr, 3*ct.num_ions, MPI_DOUBLE, 0, pct.img_comm );
		MPI_Bcast( R_ptr, 3*ct.num_ions, MPI_DOUBLE, 0, pct.img_comm );
		MPI_Bcast( &L_total, 1, MPI_DOUBLE, 0, pct.img_comm );
		MPI_Bcast( &R_total, 1, MPI_DOUBLE, 0, pct.img_comm );

        //dprintf(" %d: Do(n't) apply constraints to first/last images ", constrain);
        MPI_Bcast( &constrain, 1, MPI_INT, 0, pct.img_comm );
        ct.constrainforces = constrain;

		/* capture force constraint parameters from right and left data*/
		for ( count = 0; count < ct.num_ions; count++ )
		{
			/* put force constraints into control structure */
			ct.ions[count].constraint.setA_weight = L_total;
			ct.ions[count].constraint.setA_coord[X] = L_ptr[3*count + X]; 
			ct.ions[count].constraint.setA_coord[Y] = L_ptr[3*count + Y]; 
			ct.ions[count].constraint.setA_coord[Z] = L_ptr[3*count + Z]; 

			ct.ions[count].constraint.setB_weight = R_total;
			ct.ions[count].constraint.setB_coord[X] = R_ptr[3*count + X]; 
			ct.ions[count].constraint.setB_coord[Y] = R_ptr[3*count + Y]; 
			ct.ions[count].constraint.setB_coord[Z] = R_ptr[3*count + Z]; 

			/* zero velocities for every nudge */
			ct.ions[count].velocity[0] = 0.0;
			ct.ions[count].velocity[1] = 0.0;
			ct.ions[count].velocity[2] = 0.0;
		}

        /* Call fastrelax for max_md_steps steps */
		MPI_Barrier( MPI_COMM_WORLD );
		printf("\tNEB call fast relax.\n");
		fflush(NULL);fsync( fileno(ct.logfile) );
		MPI_Allreduce( &tmp_mag, &max_frc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        relax (2, states, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc);

        /* Check for NEB convergence */
        /* Are we force converged? */
        tmp_mag = 0.0;
        max_frc = 0.0;
        for (count = 0; count < ct.num_ions; count++)
        {
            fp = ct.ions[count].force[ct.fpt[0]];
            tmp_mag =  fp[X] * fp[X] + fp[Y] * fp[Y] + fp[Z] * fp[Z];
            if ( tmp_mag > max_frc )
                max_frc = tmp_mag;
        }

        printf(" Find max force amongst all images ");
		fflush(NULL);fsync( fileno(ct.logfile) );

		tmp_mag = max_frc;

		MPI_Barrier( MPI_COMM_WORLD );
        printf(" Passed Global max force barrier ");

		MPI_Allreduce( &tmp_mag, &max_frc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        printf(" Max force to all procs, loop termination condition ");
		
        MPI_Bcast( &max_frc, 1, MPI_DOUBLE, 0, pct.img_comm );

        CONV_FORCE = (max_frc < ct.thr_frc * ct.thr_frc);
		printf("\nNEB is max_frc^2:%f < ct.thr_frc^2:%f ? If so, DONE.\n", max_frc, ct.thr_frc*ct.thr_frc);
		printf("\nNEB is neb_steps:%d== max_rmg_steps:%d ? If so, DONE.\n", neb_steps, ct.max_rmg_steps);

        DONE = (CONV_FORCE || (++neb_steps == ct.max_rmg_steps));

    } /* end while(!DONE) */
}                               /* end neb_relax */


/******/
