/************************** SVN Revision Information **************************
 **    $Id: get_te1.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*

                        get_te.c


    Gets total energy of the system. 
    Stores result in control structure.


*/



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"


void get_te1(double *rho, double *rhoc, double *rhocore, double *vh, double *vxc, STATE * states)
{
    double eigsum, energy_sc;
    double time1, time2;


    REAL esum, Evxcold_rho, Evhold_rho, Evh_rho, Evh_rhoc;
    REAL *vh_tem;
    int ione = 1;


    vh_tem = work_memory;
    time1 = my_crtc();


    eigsum = get_sum_eig(states);

    ct.II = get_te_ion_ion();

    ct.XC = get_Exc(rho, rhocore);

    esum = ddot(&P0_BASIS, vxc, &ione, rho, &ione);
    Evxcold_rho = ct.vel * real_sum_all(esum);
    esum = ddot(&P0_BASIS, vh, &ione, rho, &ione);
    Evhold_rho = 0.5 * ct.vel * real_sum_all(esum);

    get_vh(rho, rhoc, vh_tem, 10, ct.poi_parm.levels);
    esum = ddot(&P0_BASIS, vh_tem, &ione, rho, &ione);
    Evh_rho = 0.5 * ct.vel * real_sum_all(esum);
    esum = ddot(&P0_BASIS, vh_tem, &ione, rhoc, &ione);
    Evh_rhoc = 0.5 * ct.vel * real_sum_all(esum);


    /* SC energy correction: the eigenvalues are computed with
       the "old" potential, ct.XC with the the new one */
    energy_sc = eigsum + ct.XC + ct.II - Evxcold_rho - 2. * Evhold_rho + Evh_rho - Evh_rhoc;
    ct.TOTAL = energy_sc;

    if (pct.thispe == 0)
    {

        printf("\n\n\n @@ EIGENVALUE SUM = %13.5f", eigsum);
        printf("\n @@ ION_ION            = %13.7f", ct.II);
        printf("\n @@ ELECTROSTATIC      = %13.7f", Evh_rho);
        printf("\n @@ ELECTROSTATIC rhoc = %13.7f", Evh_rhoc);
        printf("\n @@ ELECTROSTATIC old  = %13.7f", Evhold_rho);
        printf("\n @@ XC                 = %13.7f", ct.XC);
        printf("\n @@ XC potential       = %13.7f", Evxcold_rho);
        printf("\n @@ SC Energy  aaa     = %13.7f", ct.TOTAL);

    }                           /* end if */



    time2 = my_crtc();
    rmg_timings(GET_TE_TIME, (time2 - time1), 0);


}                               /* end get_te */
