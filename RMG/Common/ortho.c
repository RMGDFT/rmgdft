/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <stdio.h>
#include <float.h>
#include <math.h>


#if FAST_ORTHO
/*Orthogonalizes states for ALL k-points
 * ortho is a similar function, except that it does orthogonalization for specified k-point*/
void ortho (STATE * states, int kpt)
{
    int ist1, ist2, length, size, ione=1, idx;
    REAL time1, time2, rone=1.0;
    REAL *cR, *cI, *Oij;
    STATE *st, *st1;

    if(ct.norm_conserving_pp) {
        ortho_ncpp (states);
        return;
    }

    time1 = my_crtc ();

    my_malloc (cR, ct.num_states, REAL);
    my_malloc (cI, ct.num_states, REAL);
    my_malloc (Oij, ct.num_states * ct.num_states, REAL);


    st = states;

#if MD_TIMERS
    time2 = my_crtc ();
#endif
    for (ist1 = 0; ist1 < ct.num_states; ist1++) {
        st1 = &st[ist1];
        norm_psi1 (st1, ist1, kpt);
    }
#if MD_TIMERS
    rmg_timings (ORTHO_NORM_PSI, (my_crtc () - time2));
    time2 = my_crtc ();
#endif

    size =pct.P0_BASIS;
    get_psi_overlaps(st->psiR, Oij, ct.num_states, ct.num_states, size, size);
#if MD_TIMERS
    rmg_timings (ORTHO_GET_OVERLAPS, (my_crtc () - time2));
    time2 = my_crtc ();
#endif


    for (ist1 = 0; ist1 < ct.num_states; ist1++)
    {

        st1 = &st[ist1];

#if MD_TIMERS
        time2 = my_crtc ();
#endif
        /*This will calculate cR and cI coefficients */
        for (ist2 = ist1 + 1; ist2 < ct.num_states; ist2++)
            ortho_get_coeff (st1, &st[ist2], ist1, ist2, kpt, &cR[ist2], &cI[ist2], Oij);
#if MD_TIMERS
        rmg_timings (ORTHO_GET_COEFF, (my_crtc () - time2));
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
        rmg_timings (ORTHO_GLOB_SUM, (my_crtc () - time2));
        time2 = my_crtc ();
#endif

#if 0
        /*Update wavefunctions */
        for (ist2 = ist1 + 1; ist2 < ct.num_states; ist2++) {
            /*update the wavefunction psi2 */
            QMD_daxpy (size, cR[ist2], st1->psiR, ione, st[ist2].psiR, ione);
        }
#endif
        idx = ct.num_states - ist1 - 1;
        if(idx)
            dger_(&size, &idx, &rone, st[ist1].psiR, &ione,
                    &cR[ist1+1], &ione, st[ist1+1].psiR, &size);


        for (ist2 = ist1 + 1; ist2 < ct.num_states; ist2++) {
            update_waves (st1, &st[ist2], ist1, ist2, kpt, cR[ist2], cI[ist2]);

        }
        norm_psi1 (st1, ist1, kpt);

#if MD_TIMERS
        rmg_timings (ORTHO_UPDATE_WAVES, (my_crtc () - time2));
#endif
    }                       /*end for ist1 */



    my_free (Oij);
    my_free (cI);
    my_free (cR);

    rmg_timings (ORTHO_TIME, (my_crtc () - time1));

}                               /*end ortho_full */

#else


/*Orthogonalizes states for ALL k-points
 * ortho is a similar function, except that it does orthogonalization for specified k-point*/
void ortho (STATE * states, int kpt)
{
    int ist1, ist2, length;
    REAL time1, time2;
    REAL *cR, *cI;
    STATE *st, *st1;

    if(ct.norm_conserving_pp) {
        ortho_ncpp (states);
        return;
    }

    time1 = my_crtc ();

    my_malloc (cR, ct.num_states, REAL);
    my_malloc (cI, ct.num_states, REAL);




    st = states;
    for (ist1 = 0; ist1 < ct.num_states; ist1++)
    {

        st1 = &st[ist1];

#if MD_TIMERS
        time2 = my_crtc ();
#endif
        norm_psi1 (st1, ist1, kpt);
#if MD_TIMERS
        rmg_timings (ORTHO_NORM_PSI, (my_crtc () - time2));
        time2 = my_crtc ();
#endif

        /*This will calculate cR and cI coefficients */
        for (ist2 = ist1 + 1; ist2 < ct.num_states; ist2++)
            ortho_get_coeff (st1, &st[ist2], ist1, ist2, kpt, &cR[ist2], &cI[ist2]);

#if MD_TIMERS
        rmg_timings (ORTHO_GET_COEFF, (my_crtc () - time2));
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
        rmg_timings (ORTHO_GLOB_SUM, (my_crtc () - time2));
        time2 = my_crtc ();
#endif

        /*Update wavefunctions */
        for (ist2 = ist1 + 1; ist2 < ct.num_states; ist2++)
            update_waves (st1, &st[ist2], ist1, ist2, kpt, cR[ist2], cI[ist2]);

#if MD_TIMERS
        rmg_timings (ORTHO_UPDATE_WAVES, (my_crtc () - time2));
#endif
    }                       /*end for ist1 */


    my_free (cR);
    my_free (cI);

    rmg_timings (ORTHO_TIME, (my_crtc () - time1));

}                               /*end ortho_full */
#endif
