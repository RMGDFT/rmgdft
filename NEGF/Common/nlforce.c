/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*nlforce.c *****
 */



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"



void nlforce (REAL * veff)
{
    int ion, i, isp, count;
    int nh, size, n2, idx;
    REAL *rho_nm, *QnmI_R, *forces_tem;
    REAL *part_rho_nm_x, *part_rho_nm_y, *part_rho_nm_z;
    ION *iptr;
    SPECIES *sp;

    REAL time1, time2;
    time1 = my_crtc();

 
    size = ct.num_ions * ct.max_nl * ct.max_nl;
    my_malloc_init( rho_nm, 4 * size, REAL );
    my_malloc_init( forces_tem, ct.num_ions*3, REAL );
    part_rho_nm_x = rho_nm + size;
    part_rho_nm_y = part_rho_nm_x + size;
    part_rho_nm_z = part_rho_nm_y + size;

    get_all_partial_kbpsi (states);

    n2 = ct.num_states * ct.num_states;

    tri_to_row (lcr[0].density_matrix_tri, work_matrix, ct.num_blocks, ct.block_dim);
    rho_Qnm_mat (rho_nm, work_matrix);
    partial_Mat_nm_R (part_rho_nm_x, part_rho_nm_y, part_rho_nm_z, work_matrix);

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];
        nh = sp->num_projectors;

        nlforce_par_Q(veff, rho_nm, ion, nh, &forces_tem[ion*3]);

        /*
           if(pct.gridpe == 0) {
           printf("\nafter nlforce_par_Q \n");
           write_force();
        }

        nlforce_par_rho(part_rho_nm_x, part_rho_nm_y, part_rho_nm_z, ion, nh);
        /*
           if(pct.gridpe == 0) {
           printf("after nlforce_par_rho \n");
           write_force();
           }
         */
      
      
       
    } /* end for ion loop */

    size = 3 * ct.num_ions;
    global_sums(forces_tem, &size, pct.grid_comm);

   
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        iptr->force[ct.fpt[0]][0] += ct.vel_f * forces_tem[ion*3];
        iptr->force[ct.fpt[0]][1] += ct.vel_f * forces_tem[ion*3+1];
        iptr->force[ct.fpt[0]][2] += ct.vel_f * forces_tem[ion*3+2];
    }

    /* Initialize Non-local operators */
    is_vloc_state_overlap (states);

    init_loc_xyz ();
    get_ion_orbit_overlap_loc (states);
    partial_vloc ();

    nlforce_par_D (states, forces_tem);

    size = ct.num_ions * 3;
    global_sums (forces_tem, &size, pct.grid_comm);

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        iptr->force[ct.fpt[0]][0] += forces_tem[3*ion];
        iptr->force[ct.fpt[0]][1] += forces_tem[3*ion+1];
        iptr->force[ct.fpt[0]][2] += forces_tem[3*ion+2];
    }
	my_free (rho_nm);
    my_free (forces_tem);

    time2 = my_crtc();
    rmg_timings(NLFORCE_TIME, time2 - time1);

}
