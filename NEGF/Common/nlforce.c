/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*nlforce.c *****
 */



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"



void nlforce (REAL * veff)
{
    int ion, i, isp, count;
    int nh, size, n2, idx;
    REAL *rho_nm, *QnmI_R, *force;
    REAL *part_rho_nm_x, *part_rho_nm_y, *part_rho_nm_z;
    ION *iptr;
    SPECIES *sp;

    REAL time1, time2;
    time1 = my_crtc();

    my_malloc_init( force, 3 * ct.num_ions, REAL);

    size = ct.num_ions * ct.max_nl * ct.max_nl;
    my_malloc_init( rho_nm, 4 * size, REAL );
    part_rho_nm_x = rho_nm + size;
    part_rho_nm_y = part_rho_nm_x + size;
    part_rho_nm_z = part_rho_nm_y + size;

    get_all_partial_kbpsi (states);

    n2 = ct.num_states * ct.num_states;

    tri_to_whole_real (lcr[0].density_matrix_tri, work_matrix, ct.num_blocks, ct.block_dim);
    rho_nm_mat (rho_nm, work_matrix);
    partial_Mat_nm_R (part_rho_nm_x, part_rho_nm_y, part_rho_nm_z, work_matrix);

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];
        nh = sp->num_projectors;

        count = pct.Qidxptrlen[ion];

        if (count)
        {
            size = 3 * (nh * (nh + 1) / 2) * count;
            my_malloc_init( QnmI_R, size, REAL );

            partial_QI(ion, QnmI_R);
        }

        get_partial_ddd (QnmI_R, veff, ion, nh);

        if (count) my_free(QnmI_R);

        /*
        ********** 1. calculate -int{dr Veff(r) * sum {partial_Qnm^I(r) * rho_nm^I} }*************
        */ 
        nlforce_par_Q (rho_nm, ion, nh, force);

        /*
        ********** 2. calculate -sum {Dnm^I * partial_rho_nm^I} *************************
        */ 
        nlforce_par_rho (part_rho_nm_x, part_rho_nm_y, part_rho_nm_z, ion, nh, force);

    } /* end for ion loop */


    my_free (rho_nm);

    /*
    ************** 3. calculate -sum_{munu}( H_{munu} * partial_D_{munu}/partial_R )**********
    */

    /* Initialize Non-local operators */
    is_vloc_state_overlap (states);

    init_loc_xyz ();
    get_ion_orbit_overlap_loc (states);
    partial_vloc ();

    nlforce_par_D (states, force);

    size = ct.num_ions * 3;
    global_sums (force, &size, pct.grid_comm);

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        iptr->force[ct.fpt[0]][0] += force[3*ion];
        iptr->force[ct.fpt[0]][1] += force[3*ion+1];
        iptr->force[ct.fpt[0]][2] += force[3*ion+2];
    }

    my_free (force);

    time2 = my_crtc();
    rmg_timings(NLFORCE_TIME, time2 - time1);

}
