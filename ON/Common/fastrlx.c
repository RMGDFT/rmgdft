/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/fastrlx.c *****
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
 *   void fastrlx(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
 *                REAL *rho, REAL *rhocore, REAL *rhoc)
 *   drive routine for fast relax.
 * INPUTS
 *   states: all wave functions (see md.h)
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
 *   md.c
 * CHILDREN
 *   quench.c to_crystal.c to_cartesian.c init_nuc.c get_nlop.c scf.c
 *   get_te.c force.c md_fastrelax.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "md.h"



/* Local function prototypes */
void md_fastrelax(void);
void movie(FILE *);



void fastrlx(STATE * states, STATE * states1, REAL * vxc, REAL * vh, REAL * vnuc, REAL * vh_old,
             REAL * vxc_old, REAL * rho, REAL * rhocore, REAL * rhoc)
{

    int it, ion, level;
    int CONV_FORCE;
    int DONE;
    FILE *mfp = NULL;
    FILE *xbsfp1 = NULL;
    char xbs_filename[60];

    /* if ( ct.override_atoms == 1 ) quench(states, vxc, vh, vnuc, rho, rhocore, rhoc); */

    /* open movie file and output initial frame */
    if ((ct.rmvmovie) && (ct.max_md_steps > 0 && pct.gridpe == 0))
    {
        mfp = fopen("traj.rmv", "w");
        if (setvbuf(mfp, (char *) NULL, _IOFBF, 4096 * 16) != 0)
            printf("\n Warning: cant allocate movie io buffer size\n");

        movie(mfp);
    }

    /* open XBS movie file */
    if ((ct.xbsmovie) && (ct.max_md_steps > 0 && pct.gridpe == 0))
    {

        strcpy(xbs_filename, "traj");
        xbsfp1 = (FILE *) open_xbs_movie(xbs_filename);
    }


    /* reset total number of scf steps */
    ct.total_scf_steps = 0;


    /* ---------- begin relax loop --------- */

    /* reset number of md steps */
    ct.md_steps = 0;
    DONE = FALSE;
    while (!DONE)
    {

        if (pct.gridpe == 0)
            printf("\nfastrlx: ---------- [md: %d/%d] ----------\n", ct.md_steps, ct.max_md_steps);

        /* quench the electrons and calculate forces */
        quench(states, states1, vxc, vh, vnuc, vh_old, vxc_old, rho, rhoc, rhocore);


        /* to_crystal enforces periodic boundary conditions */
        for (it = 0; it < ct.num_ions; it++)
        {
            to_crystal(ct.ions[it].xtal, ct.ions[it].crds);
            to_cartesian(ct.ions[it].xtal, ct.ions[it].crds);
        }

        /* output the milliken charges */
/*    if (ct.domilliken == 1) {
      ct.milliken = 1;
      get_milliken(states);
      get_nlop();
    }
*/

        /* save data to file for future restart */
        if (ct.checkpoint)
            if ((ct.md_steps % ct.checkpoint == 0) && (ct.md_steps))
                write_data(ct.outfile, vh, vxc, vh_old, vxc_old, rho, states);


        /* check force convergence */
        CONV_FORCE = TRUE;
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            if (ct.ions[ion].movable)
            {
                REAL *fp;
                fp = ct.ions[ion].force[ct.fpt[0]];
                CONV_FORCE &= ((fabs(fp[0]) < ct.thr_frc) &&
                               (fabs(fp[1]) < ct.thr_frc) && (fabs(fp[2]) < ct.thr_frc));
            }
        }

        /* done if forces converged or reached limit of md steps */
        DONE = (CONV_FORCE || ct.md_steps == ct.max_md_steps);





        /* not done yet ? => move atoms */
        if (!DONE)
        {

            /* move the ions and update the coordinates */
            md_fastrelax();

            /* update the coordinates of the states corresponding the ions */
            change_states_crds(states);

            state_corner_xyz(states);

            /* determine the overlap regions */
            is_state_overlap(states);
            get_orbit_overlap_region(states);

            /* set up the communication matrix */
            init_comm(states);
            init_comm_res(states);

            /* duplicate the information of states */
            duplicate_states_info(states, states1);
            duplicate_states_info(states, states_tem);

            for (level = 0; level < ct.eig_parm.levels + 1; level++)
                make_mask_grid_state(level, states);

            /* Initialize the nuclear local potential and the compensating charges */
            init_nuc(vnuc, rhoc, rhocore);

            /* Initialize Non-local operators */
            init_nl_xyz();
            get_ion_orbit_overlap_nl(states);
            get_nlop();
            init_nonlocal_comm();

            /* Initialize qfuction in Cartesin coordinates */
            init_qfunct();
            get_QI();

            /* Get the qqq */
            get_qqq();


            /* write out frame */
            if ((ct.rmvmovie) && (ct.max_md_steps > 0 && pct.gridpe == 0))
                movie(mfp);


            /* output xbs frame */
            if ((ct.xbsmovie) && (ct.max_md_steps > 0 && pct.gridpe == 0))
            {

                if (ct.md_steps % ct.xbsmovie == 0)
                    xbsmovie(xbsfp1);

                /*Flush the file from time to time */
                if ((ct.md_steps) && (ct.md_steps % 10 == 0))
                    fflush(xbsfp1);
            }                   /*end if ((ct.xbsmovie) && (ct.max_md_steps > 0 && pct.gridpe == 0)) */

            /* ct.md_steps measures the number of updates to the atomic positions */
            ct.md_steps++;

        }                       /*end if (!DONE) */


    }
    /* ---------- end relax loop --------- */

    if (ct.max_md_steps > 0 && pct.gridpe == 0)
    {

        printf("\n");
        progress_tag();

        if (CONV_FORCE)
            printf("force convergence has been achieved. stopping ...\n");
        else
            printf
                ("force convergence has NOT been achieved. stopping (max number of MD steps reached) ...\n");

    }



    /*Write out final data */
    /*write_data(ct.outfile, vh, vxc, vh_old, vxc_old, rho, states); */

    /* close moviefile */
    if ((ct.rmvmovie) && (ct.max_md_steps > 0 && pct.gridpe == 0))
        if ((ct.xbsmovie) && (ct.max_md_steps > 0 && pct.gridpe == 0))
            fclose(mfp);


}                               /* end fastrlx */


/******/
