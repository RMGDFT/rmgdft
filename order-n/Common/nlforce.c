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



void nlforce(REAL * veff)
{
    int ion, i, isp;
    int nh, size, n2, idx;
    REAL *rho_nm;
    REAL *part_rho_nm_x, *part_omega_nm_x;
    REAL *part_rho_nm_y, *part_omega_nm_y;
    REAL *part_rho_nm_z, *part_omega_nm_z;
    ION *iptr;
    SPECIES *sp;
    double *forces_tem;

    int IA=1,JA=1,IB=1,JB=1, numst = ct.num_states;

    REAL time1, time2;
    time1 = my_crtc();

    size = ct.num_ions * ct.max_nl * ct.max_nl;
    my_malloc_init( rho_nm, 7 * size, REAL );
    my_malloc_init( forces_tem, ct.num_ions*3, REAL );
    part_rho_nm_x = rho_nm + size;
    part_rho_nm_y = part_rho_nm_x + size;
    part_rho_nm_z = part_rho_nm_y + size;
    part_omega_nm_x = part_rho_nm_z + size;
    part_omega_nm_y = part_omega_nm_x + size;
    part_omega_nm_z = part_omega_nm_y + size;

    for (idx = 0; idx < 7 * size; idx++)
    {
        rho_nm[idx] = 0.0;
    }

    get_all_partial_kbpsi(states);

    n2 = ct.num_states * ct.num_states;

    Cpdgemr2d(numst, numst, mat_X, IA, JA, pct.desca, work_matrix_row,
            IB, JB, pct.descb, pct.descb[1]);


    rho_nm_mat(rho_nm, work_matrix_row);
    partial_Mat_nm_R(part_rho_nm_x, part_rho_nm_y, part_rho_nm_z, work_matrix_row);

    Cpdgemr2d(numst, numst, mat_Omega, IA, JA, pct.desca, work_matrix_row,
            IB, JB, pct.descb, pct.descb[1]);
    partial_Mat_nm_R(part_omega_nm_x, part_omega_nm_y, part_omega_nm_z, work_matrix_row);

    global_sums(rho_nm, &size, pct.grid_comm);
    global_sums(part_rho_nm_x, &size, pct.grid_comm);
    global_sums(part_rho_nm_y, &size, pct.grid_comm);
    global_sums(part_rho_nm_z, &size, pct.grid_comm);
    global_sums(part_omega_nm_x, &size, pct.grid_comm);
    global_sums(part_omega_nm_y, &size, pct.grid_comm);
    global_sums(part_omega_nm_z, &size, pct.grid_comm);


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
         */

        nlforce_par_rho(part_rho_nm_x, part_rho_nm_y, part_rho_nm_z, ion, nh);
        /*
           if(pct.gridpe == 0) {
           printf("after nlforce_par_rho \n");
           write_force();
           }
         */

        nlforce_par_omega(part_omega_nm_x, part_omega_nm_y, part_omega_nm_z, ion, nh);
        /*
           if(pct.gridpe == 0) {
           printf("\nafter nlforce_par_omega \n");
           write_force();
           }
         */

    }

    size = 3 * ct.num_ions;
    global_sums(forces_tem, &size, pct.grid_comm);

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        iptr->force[ct.fpt[0]][0] += ct.vel_f * forces_tem[ion*3];
        iptr->force[ct.fpt[0]][1] += ct.vel_f * forces_tem[ion*3+1];
        iptr->force[ct.fpt[0]][2] += ct.vel_f * forces_tem[ion*3+2];
    }




    my_free(rho_nm);
    my_free(forces_tem);

    time2 = my_crtc();
    rmg_timings(NLFORCE_TIME, time2 - time1);

}


