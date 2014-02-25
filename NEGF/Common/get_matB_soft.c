/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"

/* 
   output: matB = matrix (psi**T * B * psi) 
*/






#if USE_DIS_MAT

#include "my_scalapack.h"

#endif /* USE_DIS_MAT */


void get_matB_soft (STATE * states, STATE * states1, double *mat)
{
    int istate, jstate;
    int n2 = ct.num_states * ct.num_states;
    int ione = 1;
    int maxst = ct.num_states, dim;
    double zero = 0., one = 1.;
    char trans;
    STATE state1;
    int ion, ixx, iyy, izz;
    double temp, temp1;
    double time1, time2, time3, time4;

    time1 = my_crtc ();

    for (istate = 0; istate < n2; istate++)
        mat[istate] = 0.;


    /* calculate the < states.psiR | states.psiR>  */

    my_barrier ();


    time3 = my_crtc ();

    orbit_dot_orbit (states, states, mat);

    my_barrier ();


    time4 = my_crtc ();
    rmg_timings (ORBIT_DOT_ORBIT_O, (time4 - time3));


        /** Add sum_n,m,I(q_n,m * <phi|beta_n,I> * <beta_m,I|phi>) **/
    get_matB_qnm (mat);         /* shuchun wang */

    time3 = my_crtc ();
    rmg_timings (matB_qnm_TIME, (time3 - time4));


    /* Sum over all processors */
    global_sums (mat, &n2, pct.grid_comm);


    dscal (&n2, &ct.vel, mat, &ione);

    /* symmetrize the mat */

    for (istate = 0; istate < ct.num_states - 1; istate++)
    {
        for (jstate = istate + 1; jstate < ct.num_states; jstate++)
        {
            mat[istate * ct.num_states + jstate] = 0.5 * (mat[istate * ct.num_states + jstate] +
                                                          mat[jstate * ct.num_states + istate]);
            mat[jstate * ct.num_states + istate] = mat[istate * ct.num_states + jstate];
        }
    }

    if (pct.gridpe == 0)
    {
        print_matrix (mat, 5, maxst);
    }

    time2 = my_crtc ();
    rmg_timings (GET_MATB_SOFT_TIME, (time2 - time1));


}
