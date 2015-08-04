/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"

#include "my_scalapack.h"


void get_new_rho(STATE * states, double *rho)
{
    int ii, idx, ione = 1;
    rmg_double_t t2;
    register double tcharge;

    /* for parallel libraries */

    rmg_double_t *psi1, *psi2, *psi3, *psi_p, scale;
    int i, st1, st2, proc1, proc2;
    int loop, state_per_proc, num_send, num_recv, num_sendrecv, size1, size2;
    MPI_Status mstatus;
    rmg_double_t *rho_temp;
    MPI_Request mr_send, *mr_recv;

    int IA=1,JA=1,IB=1,JB=1, numst = ct.num_states;
    int st11;

    void *RT0 = BeginRmgTimer("3-get_new_rho");
    state_per_proc = ct.state_per_proc + 2;
    if (pct.gridpe == 0)
        printf(" Compute new density\n");

    my_malloc_init( rho_temp, get_P0_BASIS(), rmg_double_t );

#if  	DEBUG
    print_sum_square(get_P0_BASIS(), rho, "rho_sum_square before get_new_rho  ");
#endif


    my_barrier();
    void *RT = BeginRmgTimer("3-get_new_rho: Cpdgemr2d");

   Cpdgemr2d(numst, numst, mat_X, IA, JA, pct.desca, work_matrix_row, IB, JB,
            pct.descb, pct.desca[1]);

    EndRmgTimer(RT);

   for (idx = 0; idx < get_NX_GRID() * get_NY_GRID() * get_NZ_GRID(); idx++)
       rho_global[idx] = 0.;

    void *RT1 = BeginRmgTimer("3-get_new_rho: states in this proc");
   for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
       for (st2 = st1; st2 < ct.state_end; st2++)
       {
           st11 = st1 - ct.state_begin;
           if (st1 == st2)
               scale =  work_matrix_row[st11 * ct.num_states + st2];
           if (st1 != st2) scale = 2.0 * work_matrix_row[st11 * ct.num_states + st2];
           psi1 = states[st1].psiR;
           psi2 = states[st2].psiR;

           if (state_overlap_or_not[st11 * ct.num_states + st2 ] == 1)
               density_orbit_X_orbit(st1, st2, scale, psi1, psi2,
                       rho_global, 0, states, orbit_overlap_region);

       }

    EndRmgTimer(RT1);

    void *RT2 = BeginRmgTimer("3-get_new_rho: states other proc");
   my_malloc_init(psi3, ct.max_orbit_size, rmg_double_t );

   for (loop = 0; loop < num_sendrecv_loop; loop++)
   {
       proc2 = recv_from[loop * state_per_proc];
       num_recv = recv_from[loop * state_per_proc + 1];

       for (i = 0; i < num_recv; i++)
       {
           st2 = recv_from[loop * state_per_proc + i + 2];

           psi2= states[st2].psiR;
           for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
           {
               psi1 = states[st1].psiR;

               st11 = st1 - ct.state_begin;

               if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
               {
                   psi1 = states[st1].psiR;
                   scale = 2.0 * work_matrix_row[st11 * ct.num_states + st2];
                   density_orbit_X_orbit(st1, st2, scale, psi1, psi2,
                           rho_global, 0, states, orbit_overlap_region);
               }
           }
       }

   }                           /* end of loop  */

   my_barrier();

    EndRmgTimer(RT2);

    void *RT3 = BeginRmgTimer("3-get_new_rho: distribution");
   idx = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
   global_sums(rho_global, &idx, pct.grid_comm);

   global_to_distribute(rho_global, rho_temp);

    EndRmgTimer(RT3);

    void *RT4 = BeginRmgTimer("3-get_new_rho: interpolation");
   mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);

   my_free(rho_temp);

    EndRmgTimer(RT4);

    void *RT5 = BeginRmgTimer("3-get_new_rho: augmented");

   rho_augmented(rho, work_matrix_row, state_begin, state_end, num_nonlocal_ion, 
           kbpsi, max_ion_nonlocal, kbpsi_comm, ionidx_allproc);

    EndRmgTimer(RT5);

   int iii = get_FP0_BASIS();

   tcharge = 0.0;
   for (idx = 0; idx < get_FP0_BASIS(); idx++)
       tcharge += rho[idx];
   ct.tcharge = real_sum_all(tcharge, pct.grid_comm);
   ct.tcharge = real_sum_all(ct.tcharge, pct.spin_comm);


   ct.tcharge *= get_vel_f();

   t2 = ct.nel / ct.tcharge;
   sscal(&iii, &t2, &rho[0], &ione);
    

   if (pct.gridpe == 0)
       printf("\n total charge Normalization constant = %f  \n", t2);

    EndRmgTimer(RT0);


#if  	DEBUG
   print_sum_square(get_P0_BASIS(), rho, "rho_sum_sqare in the end of get_new_rho  ");
#endif

}
