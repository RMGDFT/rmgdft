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

#include <semaphore.h>

#include "grid.h"
#include "rmgtypes.h"
#include "common_prototypes.h"
#include "main.h"
#include "rmgthreads.h"


void init_pe ( int image )
{

    int ii, jj, kk, ix, iy, iz, idx, ioffset, rem, ierr, thread;
    int image_grp_map[MAX_IMGS], range[1][3], neighbors[6];
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


    // Not sure if this will work on anything other than the Crays ....
    if(ct.images_per_node > 1) {
    
        int ix, ir1, ir2, ir3, img_nxdim, img_nydim, img_nx, img_ny, worldpe;
        MPI_Comm_rank (MPI_COMM_WORLD, &worldpe);

        img_nxdim = pct.images / ct.images_per_node;
        img_nydim = ct.images_per_node;
        img_nx = pct.thisimg / img_nydim;
        img_ny = pct.thisimg % img_nydim;
        range[0][0] = img_nx * img_nydim * pct.grids * NPES + img_ny;
        range[0][1] = img_nx * img_nydim * pct.grids * NPES + img_ny + img_nydim*NPES - img_nydim;
        range[0][2] = ct.images_per_node;
        Dprintf("DEBUG %d  %d  %d  %d  %d  %d  %d",pct.thisimg, img_nx, img_ny, pct.grids, range[0][0],range[0][1],range[0][2]);sleep(2);
        /* define this images group and put its comm in pct */
        ierr=MPI_Group_range_incl (world_grp, 1, range, &group);
        Dprintf("IERR0 = %d  WPE=%d",ierr, worldpe);fflush(NULL);
        ierr=MPI_Comm_create (MPI_COMM_WORLD, group, &pct.img_comm);
        Dprintf("IERR1 = %d  WPE=%d",ierr, worldpe);fflush(NULL);

        /* setup pct.rmg_comm to include all image group_rank 0 pe's */
        /* build rank list of group masters, this assumes contiguous ranges */
        /* NOTE: this explicitly depends on range assignment method above! */
// Still needs to be extended for images stacked horizontally
        for (image = 0; image < pct.images; image++) {
            img_nx = image / ct.images_per_node; 
            img_ny = image % ct.images_per_node;
            image_grp_map[image] = img_nx * img_nydim * pct.grids * NPES + img_ny;
        }

    }
    else {

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

    }


    MPI_Barrier(MPI_COMM_WORLD);
    /* define master group and make its comm global */
    MPI_Group_incl (world_grp, pct.images, image_grp_map, &img_masters);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_create (MPI_COMM_WORLD, img_masters, &pct.rmg_comm);
    MPI_Barrier(MPI_COMM_WORLD);

    /* set gridpe rank value to its image rank */
    ierr = MPI_Comm_rank (pct.img_comm, &pct.imgpe);

    /* Read in our control information, depends on pct.img_comm for dissemination */
    read_control (ct.cfile);

    /* this will need to be extended if we want to parallelize over k_points */
    if (pct.grids == 1)
    {
        if (ct.spin_flag)
            error_handler
                ("Spin calculations require 2 grids, please rerun with twice as many PEs.");
        MPI_Barrier(MPI_COMM_WORLD);
        ierr = MPI_Comm_dup (pct.img_comm, &pct.grid_comm);
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

    MPI_Barrier(MPI_COMM_WORLD);
    /* set gridpe rank value to local grid rank value */
    MPI_Comm_rank (pct.grid_comm, &pct.gridpe);
    Dprintf("My grid rank is %d and my image rank is %d\n", pct.gridpe, pct.imgpe);

    // Set up grids and neighbors
    set_rank(pct.gridpe);

}                               /* end init_pe */
