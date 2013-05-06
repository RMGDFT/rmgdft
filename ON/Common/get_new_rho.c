/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"

#include "my_scalapack.h"


void get_new_rho(STATE * states, double *rho)
{
    int ii, idx, ione = 1;
    REAL t2;
    register double tcharge;

    /* for parallel libraries */
    int n2 = (ct.state_end-ct.state_begin)* ct.num_states;
    int mxllda;

    REAL *psi1, *psi2, *psi3, *psi_p, scale;
    int i, st1, st2, proc1, proc2;
    REAL time1, time2, time3;
    int loop, state_per_proc, num_send, num_recv, num_sendrecv, size1, size2;
    MPI_Status mstatus;
    REAL *rho_temp;
    MPI_Request mr_send, *mr_recv;

    int IA=1,JA=1,IB=1,JB=1, numst = ct.num_states;
    int st11;

    state_per_proc = ct.state_per_proc + 2;
    time1 = my_crtc();
    if (pct.gridpe == 0)
        printf(" Compute new density\n");

    my_malloc_init( rho_temp, P0_BASIS, REAL );

#if  	DEBUG
    print_sum_square(P0_BASIS, rho, "rho_sum_square before get_new_rho  ");
#endif

    mxllda = MXLLDA;



   Cpdgemr2d(numst, numst, mat_X, IA, JA, pct.desca, work_matrix_row, IB, JB,
            pct.descb, pct.desca[1]);


   for (idx = 0; idx < NX_GRID * NY_GRID * NZ_GRID; idx++)
       rho_global[idx] = 0.;

   time2 = my_crtc();
   for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
       for (st2 = st1; st2 < ct.state_end; st2++)
       {
           st11 = st1 - ct.state_begin;
           if (st1 == st2)
               scale =  work_matrix_row[st11 * ct.num_states + st2];
           if (st1 != st2) scale = 2.0 * work_matrix_row[st11 * ct.num_states + st2];
           psi1 = states[st1].psiR;
           psi2 = states[st2].psiR;

           if (state_overlap_or_not[st1 + st2 * ct.num_states] == 1)
               density_orbit_X_orbit(st1, st2, scale, psi1, psi2, rho_global, 0, states);

       }

   int max_ii = 0;
   for (loop = 0; loop < num_sendrecv_loop; loop++)
   {
       num_send = send_to[loop * state_per_proc + 1];
       num_recv = recv_from[loop * state_per_proc + 1];
       max_ii = max(max_ii, num_send);
       max_ii = max(max_ii, num_recv);
   }

   max_ii = int_max_all(max_ii);

   ii = num_sendrecv_loop * (max_ii +10) +1;

   my_calloc(mr_recv, ii, MPI_Request);

   psi2 = orbit_tem;
   my_malloc_init(psi3, ct.max_orbit_size, REAL );

   for (loop = 0; loop < num_sendrecv_loop; loop++)
   {
       my_barrier();
       ii = loop * max_ii +1;
       proc1 = send_to[loop * state_per_proc];
       proc2 = recv_from[loop * state_per_proc];
       num_send = send_to[loop * state_per_proc + 1];
       num_recv = recv_from[loop * state_per_proc + 1];
       num_sendrecv = min(num_send, num_recv);

       if(num_sendrecv >0)
       {
           ii++;
           i=0;
           st1 = send_to[loop * state_per_proc + i + 2];
           st2 = recv_from[loop * state_per_proc + i + 2];
           size1 = states[st1].size;
           size2 = states[st2].size;
           psi1 = states[st1].psiR;

           MPI_Isend(psi1, size1, MPI_DOUBLE, proc1, ii, MPI_COMM_WORLD, &mr_send);
           MPI_Request_free(&mr_send);
           if(ii%2 == 0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
           if(ii%2 == 1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
       }



       for (i = 1; i < num_sendrecv; i++)
       {
           ii++;
           st1 = send_to[loop * state_per_proc + i + 2];
           st2 = recv_from[loop * state_per_proc + i + 2];

           psi1 = states[st1].psiR;
           size1 = states[st1].size;
           size2 = states[st2].size;

           MPI_Wait(&mr_recv[ii-1], &mstatus);

           MPI_Isend(psi1, size1, MPI_DOUBLE, proc1, ii, MPI_COMM_WORLD, &mr_send);
           MPI_Request_free(&mr_send);

           if(ii%2 == 0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
           if(ii%2 == 1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);

           st2 = recv_from[loop * state_per_proc + i-1 + 2];
           if(ii%2 == 1) psi_p = psi2;
           if(ii%2 == 0) psi_p = psi3;


           for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
               if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
               {
                   st11 = st1 - ct.state_begin;
                   psi1 = states[st1].psiR;
                   scale = 2.0 * work_matrix_row[st11 * ct.num_states + st2];
                   density_orbit_X_orbit(st1, st2, scale, psi1, psi_p, rho_global, 0, states);
               }
       }

       ii++;


       if (num_send < num_recv)
       {
           i = num_send;
           st2 = recv_from[loop * state_per_proc + i + 2];
           size2 = states[st2].size;


           if(ii%2 == 0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
           if(ii%2 == 1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
       }


       if (num_send > num_recv)
       {
           i = num_recv;
           st1 = send_to[loop * state_per_proc + i + 2];
           size1 = states[st1].size;
           psi1 = states[st1].psiR;
           MPI_Isend(psi1, size1, MPI_DOUBLE, proc1, ii, MPI_COMM_WORLD, &mr_send);
           MPI_Request_free(&mr_send);
       }

       if(num_sendrecv >0)
       {

           MPI_Wait(&mr_recv[ii-1], &mstatus);

           st2 = recv_from[loop * state_per_proc + i-1 + 2];
           if(ii%2 == 1) psi_p = psi2;
           if(ii%2 == 0) psi_p = psi3;

           for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
               if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
               {
                   st11 = st1 - ct.state_begin;
                   psi1 = states[st1].psiR;
                   scale = 2.0 * work_matrix_row[st11 * ct.num_states + st2];
                   density_orbit_X_orbit(st1, st2, scale, psi1, psi_p, rho_global, 0, states);
               }


       }



       if (num_send < num_recv)
       {
           for (i = num_send+1; i < num_recv; i++)
           {
               ii++;
               st2 = recv_from[loop * state_per_proc + i + 2];
               size2 = states[st2].size;


               MPI_Wait(&mr_recv[ii-1], &mstatus);

               if(ii%2 == 0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
               if(ii%2 == 1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);

               st2 = recv_from[loop * state_per_proc + i-1 + 2];
               if(ii%2 == 1) psi_p = psi2;
               if(ii%2 == 0) psi_p = psi3;



               for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                   if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
                   {
                       st11 = st1 - ct.state_begin;
                       psi1 = states[st1].psiR;
                       scale = 2.0 * work_matrix_row[st11 * ct.num_states + st2];
                       density_orbit_X_orbit(st1, st2, scale, psi1, psi_p, rho_global, 0, states);
                   }
           }

           ii++;
           MPI_Wait(&mr_recv[ii-1], &mstatus);
           st2 = recv_from[loop * state_per_proc + i-1 + 2];
           if(ii%2 == 1) psi_p = psi2;
           if(ii%2 == 0) psi_p = psi3;

           for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
               if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
               {
                   st11 = st1 - ct.state_begin;
                   psi1 = states[st1].psiR;
                   scale = 2.0 * work_matrix_row[st11 * ct.num_states + st2];
                   density_orbit_X_orbit(st1, st2, scale, psi1, psi_p, rho_global, 0, states);
               }


       }

       if (num_send > num_recv)
           for (i = num_recv+1; i < num_send; i++)
           {
               ii++;
               st1 = send_to[loop * state_per_proc + i + 2];
               psi1 = states[st1].psiR;
               size1 = states[st1].size;


               MPI_Isend(psi1, size1, MPI_DOUBLE, proc1, ii, MPI_COMM_WORLD, &mr_send);
               MPI_Request_free(&mr_send);

           }


   }                           /* end of loop  */

   my_free(psi3);
   my_free(mr_recv);
   my_barrier();
   time3 = my_crtc();
   rmg_timings(RHO_PHI_TIME, time3 - time2);

   idx = NX_GRID * NY_GRID * NZ_GRID;
   global_sums(rho_global, &idx, pct.grid_comm);
   time2 = my_crtc();
   rmg_timings(RHO_SUM_TIME, time2 - time3);

   global_to_distribute(rho_global, rho_temp);

   mg_prolong_MAX10 (rho, rho_temp, FPX0_GRID, FPY0_GRID, FPZ0_GRID, PX0_GRID, PY0_GRID, PZ0_GRID, FG_NX, 6);

   my_free(rho_temp);

   time3 = my_crtc();
   rmg_timings(RHO_CTOF_TIME, time3 - time2);

   rho_augmented(rho, work_matrix_row);

   time2 = my_crtc();
   rmg_timings(RHO_AUG_TIME, time2 - time3);

   tcharge = 0.0;
   for (idx = 0; idx < FP0_BASIS; idx++)
       tcharge += rho[idx];
   ct.tcharge = real_sum_all(tcharge, pct.grid_comm);


   ct.tcharge *= ct.vel_f;
   t2 = ct.nel / ct.tcharge;
   sscal(&FP0_BASIS, &t2, &rho[0], &ione);

   if (pct.gridpe == 0)
       printf("\n total charge Normalization constant = %f  \n", t2);


   time1 = my_crtc() - time1;
   rmg_timings(GET_NEW_RHO, time1);

#if  	DEBUG
   print_sum_square(P0_BASIS, rho, "rho_sum_sqare in the end of get_new_rho  ");
#endif

}
