/************************** SVN Revision Information **************************
 **    $Id: get_te.c 587 2006-08-18 22:57:56Z miro $    **
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


void get_te(double *rho, double *rhoc, double *rhocore, double *vh, double *vxc, STATE * states)
{
    double eigsum, energy_sc;
    double time1, time2;


    time1 = my_crtc();


    eigsum = get_sum_eig(states);

    ct.II = get_te_ion_ion();

    ct.XC = get_Exc(rho, rhocore);


    /* SC energy correction: the eigenvalues are computed with
       the "old" potential, ct.XC with the the new one */
    energy_sc = eigsum + ct.XC + ct.II
        - ct.Evxcold_rho - 2. * ct.Evhold_rho + ct.Evh_rho - ct.Evh_rhoc;
    ct.TOTAL = energy_sc;


    if (pct.thispe == 0)
    {

        printf("\n\n @@ EIGENVALUE SUM = %20.10f\n", eigsum);
        printf("@@ ION_ION            = %20.10f\n", ct.II);
        printf("@@ ELECTROSTATIC      = %20.10f\n", ct.Evh_rho);
        printf("@@ ELECTROSTATIC rhoc = %20.10f\n", ct.Evh_rhoc);
        printf("@@ ELECTROSTATIC old  = %20.10f\n", ct.Evhold_rho);
        printf("@@ XC                 = %20.10f\n", ct.XC);
        printf("@@ XC potential       = %20.10f\n", ct.Evxcold_rho);
        printf("@@ SC ENERGY          = %20.12f\n", energy_sc);

    }                           /* end if */


    time2 = my_crtc();
    rmg_timings(GET_TE_TIME, (time2 - time1), 0);


}                               /* end get_te */
