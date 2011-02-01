/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/****f* QMD-MGDFT/init_pe.c *****
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
 *   void init_pe(void)
 *   Initializes processor control structures.
 *   Make each processor knows who are neighbors
 * INPUTS
 *   nothing
 * OUTPUT
 *   neighbor informations are stored in pct.neighbors
 * PARENTS
 *   main.c
 * CHILDREN
 *   pe2xyz.c  
 * SOURCE
 */

#include "main.h"


void init_pe ()
{

    int npes, ii, jj, kk;
    char tmpstring[MAX_PATH], logname[MAX_PATH], basename[MAX_PATH];
    int instance, mastergrp[pct.instances], range[1][3];
    MPI_Group group, grp_world, grp_master;

    /* Setup values for this parallel run */
    MPI_Comm_size (MPI_COMM_WORLD, &npes);

    /* Check for sufficient processors for (multi-instance) run */
    if (NPES * pct.instances != npes)
        error_handler
            ("RMG compiled for %d PEs and running %d instances, but allocated %d PEs != PEs*instances !!\n",
             NPES, pct.instances, npes);

    /* setup pct.thisgrp_comm to include all pe's in this instance */
    /* get world group handle */
    MPI_Comm_group (MPI_COMM_WORLD, &grp_world);

    /* determine range of pe ranks in this instance */
    range[0][0] = pct.thisgrp * NPES;
    range[0][1] = (pct.thisgrp + 1) * NPES - 1;
    range[0][2] = 1;

    /* define this instances group and make its comm global */
    MPI_Group_range_incl (grp_world, 1, range, &group);
    MPI_Comm_create (MPI_COMM_WORLD, group, &pct.thisgrp_comm);

    /* reset thispe rank value to local instance value */
    MPI_Group_rank (group, &pct.thispe);


    /* setup pct.master_comm to include all group_rank 0 pe's */
    /* build rank list of group masters, this assumes contiguous ranges */
    /* NOTE: this explicitly depends on range assignment method above! */
	for (instance = 0; instance < pct.instances; instance++)
		mastergrp[instance] = instance * NPES;

	/* define master group and make its comm global */
	MPI_Group_incl (grp_world, pct.instances, mastergrp, &grp_master);
	MPI_Comm_create (MPI_COMM_WORLD, grp_master, &pct.master_comm);




    /* Legacy portion of init_pe */

    /* XYZ coordinates of this processor */
    pe2xyz (pct.thispe, &ii, &jj, &kk);

    /* Now wrap them in case we are running with some processors duplicated */
    /* Two should be enough for any case that we might be doing.            */
    /* Wouldn't ii %= PE_X; be better??? */
    if (ii >= PE_X)
        ii -= PE_X;
    if (ii >= PE_X)
        ii -= PE_X;

    /* Have each processor figure out who it's neighbors are */
    XYZ2PE (ii, (jj + 1) % PE_Y, kk, pct.neighbors[NB_N]);
    XYZ2PE (ii, (jj - 1 + PE_Y) % PE_Y, kk, pct.neighbors[NB_S]);
    XYZ2PE ((ii + 1) % PE_X, jj, kk, pct.neighbors[NB_E]);
    XYZ2PE ((ii - 1 + PE_X) % PE_X, jj, kk, pct.neighbors[NB_W]);
    XYZ2PE (ii, jj, (kk + 1) % PE_Z, pct.neighbors[NB_U]);
    XYZ2PE (ii, jj, (kk - 1 + PE_Z) % PE_Z, pct.neighbors[NB_D]);


    my_barrier ();
	//dprintf("Finished init_pe, all MPI groups defined.\n");

}                               /* end init_pe */
