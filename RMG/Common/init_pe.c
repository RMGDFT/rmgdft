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
 *   void init_pe( image )
 *   Initializes processor control structures. image variable
 *   contains the number of cores assigned to an image for grid setup
 *   and is reused to build the image master communicator mapping.
 *   Make each processor knows who are neighbors
 * INPUTS
 *   image
 * OUTPUT
 *   neighbor informations are stored in pct.neighbors
 * PARENTS
 *   main.c
 * CHILDREN
 *   read_control.c
 *   pe2xyz.c  
 * SOURCE
 */

#include "main.h"


void init_pe ( int image )
{

    int ii, jj, kk, ix, iy, iz, idx, ioffset, rem;
    int image_grp_map[MAX_IMGS], range[1][3];
    MPI_Group group, world_grp, img_masters;

    /* Setup MPI */
    /* get world group handle */
    MPI_Comm_group (MPI_COMM_WORLD, &world_grp);

    /* Infer the number of cpu grids in each image. */
    pct.grids = image / NPES;
    if (NPES * pct.grids != image)
        error_handler ("MPI processes per image (%d) must be a multiple of NPES (%d).", image, NPES);

    if ( pct.grids > MAX_GRIDS )
        error_handler ("CPU Grid multiplicity (%d) is more than MAX_GRIDS (%d) in params.h.", pct.grids, MAX_GRIDS);


    /* setup pct.img_comm to include all pe's in this image */

    /* determine range of pe ranks in this image */
    range[0][0] = pct.thisimg * pct.grids * NPES;
    range[0][1] = (pct.thisimg + 1) * pct.grids * NPES - 1;
    range[0][2] = 1;

    /* define this images group and put its comm in pct */
    MPI_Group_range_incl (world_grp, 1, range, &group);
    MPI_Comm_create (MPI_COMM_WORLD, group, &pct.img_comm);


    /* setup pct.rmg_comm to include all image group_rank 0 pe's */
    /* build rank list of group masters, this assumes contiguous ranges */
    /* NOTE: this explicitly depends on range assignment method above! */
    for (image = 0; image < pct.images; image++)
        image_grp_map[image] = image * NPES * pct.grids;

	/* define master group and make its comm global */
    MPI_Group_incl (world_grp, pct.images, image_grp_map, &img_masters);
    MPI_Comm_create (MPI_COMM_WORLD, img_masters, &pct.rmg_comm);

    /* set gridpe rank value to its image rank */
    MPI_Comm_rank (pct.img_comm, &pct.imgpe);

    /* Read in our control information, depends on pct.img_comm for dissemination */
    read_control (ct.cfile);

    /* this will need to be extended if we want to parallelize over k_points */
    if (pct.grids == 1)
    {
        if (ct.spin_flag)
            error_handler
                ("Spin calculations require 2 grids, please rerun with twice as many PEs.");
        MPI_Comm_dup (pct.img_comm, &pct.grid_comm);
    }
    else if (pct.grids == 2)
    {
        if (!ct.spin_flag)
            error_handler
                ("2 grids allocated but spin disabled, please rerun with half as many PEs or enable spin.");
        int ndims = 2;
        int dims[] = { NPES, 2 };
        int periods[] = { 0, 0 };
        int reorder = 1;
        int remains[] = { 1, 0 };

        MPI_Cart_create (pct.img_comm, ndims, dims, periods, reorder, &pct.grid_topo_comm);
        MPI_Cart_sub (pct.grid_topo_comm, remains, &pct.grid_comm);

        remains[0] = 0;
        remains[1] = 1;
        MPI_Cart_sub (pct.grid_topo_comm, remains, &pct.spin_comm);
    	/* set spinpe rank value to local spin rank value */
    	MPI_Comm_rank (pct.spin_comm, &pct.spinpe);

	/* for debug usage*/
	/* MPI_Comm_rank(pct.grid_comm, &pct.gridpe); */
	
	/*printf("My spin rank is %d and my image rank is %d\n", pct.spinpe, pct.imgpe);*/
	

    }
    else
    {
        error_handler ("Other than one or two grids per image not implimented.");
    }

    /* set gridpe rank value to local grid rank value */
    MPI_Comm_rank (pct.grid_comm, &pct.gridpe);
    /*printf("My grid rank is %d and my image rank is %d\n", pct.gridpe, pct.imgpe);*/

    /* Legacy portion of init_pe */

    /* XYZ coordinates of this processor */
    pe2xyz (pct.gridpe, &ii, &jj, &kk);

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

    // Compute grid sizes for each node.
    pct.PX0_GRID = NX_GRID / PE_X;
    rem = NX_GRID % PE_X;
    if(rem && (ii < rem)) pct.PX0_GRID++;

    pct.PY0_GRID = NY_GRID / PE_Y;
    rem = NY_GRID % PE_Y;
    if(rem && (jj < rem)) pct.PY0_GRID++;

    pct.PZ0_GRID = NZ_GRID / PE_Z;
    rem = NZ_GRID % PE_Z;
    if(rem && (kk < rem)) pct.PZ0_GRID++;

    pct.P0_BASIS = pct.PX0_GRID * pct.PY0_GRID * pct.PZ0_GRID;
 
    pct.FPX0_GRID = FG_NX * pct.PX0_GRID;
    pct.FPY0_GRID = FG_NY * pct.PY0_GRID;
    pct.FPZ0_GRID = FG_NZ * pct.PZ0_GRID;
    pct.FP0_BASIS = pct.FPX0_GRID * pct.FPY0_GRID * pct.FPZ0_GRID;

    // Now compute the global grid offset of the first point of the coarse and fine node grids
    find_node_offsets(pct.gridpe, NX_GRID, NY_GRID, NZ_GRID,
                      &pct.PX_OFFSET, &pct.PY_OFFSET, &pct.PZ_OFFSET);

    find_node_offsets(pct.gridpe, FNX_GRID, FNY_GRID, FNZ_GRID,
                      &pct.FPX_OFFSET, &pct.FPY_OFFSET, &pct.FPZ_OFFSET);

}                               /* end init_pe */
