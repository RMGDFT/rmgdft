/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <stdio.h>
#include <float.h>
#include <math.h>


/* Does orthogonalization for given k-point states, 
 * ortho_full is simillar it does orthoginalization for all k-points*/
void ortho (STATE * states, int kpt)
{
    int ist1, ist2, length;
    REAL time1, time2;
    REAL *cR, *cI;
    STATE *st1;

    time1 = my_crtc ();

    my_malloc (cR, ct.num_states, REAL);
    my_malloc (cI, ct.num_states, REAL);


#if MD_TIMERS
    time2 = my_crtc ();
#endif
    betaxpsi1 (states, kpt);
#if MD_TIMERS
    rmg_timings (ORTHO_BETAXPSI, (my_crtc () - time2), 0);
#endif



    for (ist1 = 0; ist1 < ct.num_states; ist1++)
    {

        st1 = &states[ist1];

#if MD_TIMERS
        time2 = my_crtc ();
#endif
        norm_psi1_parallel (st1, ist1, kpt);
#if MD_TIMERS
        rmg_timings (ORTHO_NORM_PSI, (my_crtc () - time2), 0);
        time2 = my_crtc ();
#endif

        /*This will calculate cR and cI coefficients */
        for (ist2 = ist1 + 1; ist2 < ct.num_states; ist2++)
            ortho_get_coeff (st1, &states[ist2], ist1, ist2, kpt, &cR[ist2], &cI[ist2]);

#if MD_TIMERS
        rmg_timings (ORTHO_GET_COEFF, (my_crtc () - time2), 0);
        time2 = my_crtc ();
#endif

        length = ct.num_states - (ist1 + 1);

        /*Sum coefficients over all processors */
        if (length)
        {
            global_sums (&cR[ist1 + 1], &length, pct.grid_comm);
#if !GAMMA_PT
            global_sums (&cI[ist1 + 1], &length, pct.grid_comm);
#endif
        }

#if MD_TIMERS
        rmg_timings (ORTHO_GLOB_SUM, (my_crtc () - time2), 0);
        time2 = my_crtc ();
#endif

        /*Update wavefunctions */
        for (ist2 = ist1 + 1; ist2 < ct.num_states; ist2++)
            update_waves (st1, &states[ist2], ist1, ist2, kpt, cR[ist2], cI[ist2]);

#if MD_TIMERS
        rmg_timings (ORTHO_UPDATE_WAVES, (my_crtc () - time2), 0);
#endif
    }                           /*end for ist1 */


#if MD_TIMERS
    time2 = my_crtc ();
#endif
    /*newsintR should be recalculated, since new_psi and norm_psi1_parallel do not fully
     * update newsintR so that they are efficient*/
    betaxpsi1 (states, kpt);
#if MD_TIMERS
    rmg_timings (ORTHO_BETAXPSI, (my_crtc () - time2), 0);
#endif

    my_free (cR);
    my_free (cI);

    rmg_timings (ORTHO_TIME, (my_crtc () - time1), 0);

}                               /*end ortho_full */
