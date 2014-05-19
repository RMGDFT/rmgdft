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
#include "prototypes_on.h"
#include "grid.h"


void init_pe ( int image )
{

    int ii, jj, kk;
    int image_grp_map[MAX_IMGS], range[1][3];
    MPI_Group group, world_grp, img_masters;


    /* Setup MPI */
    /* get world group handle */
    MPI_Comm_group (MPI_COMM_WORLD, &world_grp);

 /* get total mpi core count  */
    MPI_Comm_size (MPI_COMM_WORLD, &NPES);

    /* Infer the number of cpu grids in each image. */
    pct.grids = image / NPES;
    if (NPES * pct.grids != image)
        error_handler ("MPI processes per image must be a multiple of NPES (pct.pe_x*pct.pe_y*pct.pe_z).");

    if ( pct.grids > MAX_GRIDS )
        error_handler ("CPU Grid multiplicity (%d) is more than MAX_GRIDS in params.h.", pct.grids);


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



    /* this will need to be extnded if we want to parallelize over k_points */
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
    //printf("My grid rank is %d and my image rank is %d\n", pct.gridpe, pct.imgpe);

    /* Read in our control information, depends on pct.img_comm for dissemination */
    read_control (ct.cfile);




    /* Legacy portion of init_pe */

    /* XYZ coordinates of this processor */
    pe2xyz (pct.gridpe, &ii, &jj, &kk);

    /* Now wrap them in case we are running with some processors duplicated */
    /* Two should be enough for any case that we might be doing.            */
    /* Wouldn't ii %= pct.pe_x; be better??? */
    if (ii >= pct.pe_x)
        ii -= pct.pe_x;
    if (ii >= pct.pe_x)
        ii -= pct.pe_x;


    /* Have each processor figure out who it's neighbors are */
    XYZ2PE (ii, (jj + 1) % pct.pe_y, kk, pct.neighbors[NB_N]);
    XYZ2PE (ii, (jj - 1 + pct.pe_y) % pct.pe_y, kk, pct.neighbors[NB_S]);
    XYZ2PE ((ii + 1) % pct.pe_x, jj, kk, pct.neighbors[NB_E]);
    XYZ2PE ((ii - 1 + pct.pe_x) % pct.pe_x, jj, kk, pct.neighbors[NB_W]);
    XYZ2PE (ii, jj, (kk + 1) % pct.pe_z, pct.neighbors[NB_U]);
    XYZ2PE (ii, jj, (kk - 1 + pct.pe_z) % pct.pe_z, pct.neighbors[NB_D]);


 // Set up grids and neighbors using both C and C++ for now
    set_rank(pct.gridpe);


    // Compute grid sizes for each node.

//    find_node_sizes(pct.gridpe, get_NX_GRID(), get_NY_GRID(), get_NZ_GRID(), &pct.get_PX0_GRID(), &pct.get_PY0_GRID(), &pct.get_PZ0_GRID());
 //   find_node_sizes(pct.gridpe, get_FNX_GRID(), get_FNY_GRID(), get_FNZ_GRID(), &pct.get_FPX0_GRID(), &pct.get_FPY0_GRID(), &pct.get_FPZ0_GRID());

  //  pct.get_P0_BASIS() = pct.get_PX0_GRID() * pct.get_PY0_GRID() * pct.get_PZ0_GRID();
   // pct.get_FP0_BASIS() = pct.get_FPX0_GRID() * pct.get_FPY0_GRID() * pct.get_FPZ0_GRID();

    // Now compute the global grid offset of the first point of the coarse and fine node grids
    //find_node_offsets(pct.gridpe, get_NX_GRID(), get_NY_GRID(), get_NZ_GRID(),
      //                &get_PX_OFFSET(), &get_PY_OFFSET(), &get_PZ_OFFSET());

    //find_node_offsets(pct.gridpe, get_FNX_GRID(), get_FNY_GRID(), get_FNZ_GRID(),
     //                 &get_FPX_OFFSET(), &get_FPY_OFFSET(), &get_FPZ_OFFSET());


    my_barrier ();

}                               /* end init_pe */


