/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/rmg_timings.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void write_timings(void)
 *   some other subroutines to define clock for different platform
 * INPUTS
 *   nothing
 * OUTPUT
 *   print the timeings in the standard output.
 * PARENTS
 *   too many
 * CHILDREN
 *   no
 * SOURCE
 */


#include "main_on.h"
#include <time.h>
#include <stdio.h>
#include <sys/time.h>


REAL timings[LAST_TIME];


#if !HYBRID_MODEL
void rmg_timings(int what, REAL time)
{

    timings[what] += time;

}                               /* end rmg_timings */
#endif


#if (CRAY_C90)

#include <time.h>
clock_t clock(void);

REAL my_crtc(void)
{
    clock_t i;
    REAL r;

    i = clock();
    r = i * 1.0;
    return r;
}
#endif



REAL my_crtc (void)
{
    struct timeval t1;
    gettimeofday (&t1, NULL);
    return t1.tv_sec + 1e-6 * t1.tv_usec;
}



/* Outputs timing information */
void write_timings(void)
{


/* Write timing information */
    if (pct.gridpe == 0)
    {
        printf("\n\n  TIMING INFORMATION\n");
        printf("*******************************************\n");
        printf("TOTAL                          %12.4f\n", timings[TOTAL_TIME]);
        printf("-------------------------------------------\n");
        printf("  INIT                         %12.4f\n", timings[INIT_TIME]);
        printf("  QUENCH                       %12.4f\n", timings[QUENCH_TIME]);
        printf("  memory alloc.                %12.4f\n", timings[ALLOC_TIME]);
        printf("-------------------------------------------\n");
        printf("    scf                        %12.4f\n", timings[SCF_TIME]);
        printf("    force                      %12.4f\n", timings[FORCE_TIME]);
        printf("*******************************************\n");

        printf("\n");
        printf("    scf                        %12.4f\n", timings[SCF_TIME]);
        printf("-------------------------------------------\n");
        printf("      mg_eig                   %12.4f\n", timings[MG_TIME]);
        printf("      mat & diag               %12.4f\n", timings[MATDIAG_TIME]);
        printf("      get_new_rho              %12.4f\n", timings[GET_NEW_RHO]);
        printf("      update pot               %12.4f\n", timings[UPDATEPOT_TIME]);
        printf("      get_te                   %12.4f\n", timings[GET_TE_TIME]);


        printf("\n\n");
        printf("      mg_eig                   %12.4f\n", timings[MG_TIME]);
        printf("-------------------------------------------\n");
        printf("        potfc                  %12.4f\n", timings[POTFC_TIME]);
        printf("        invert matB            %12.4f\n", timings[INVMATB_TIME]);
        printf("        theta                  %12.4f\n", timings[THETA_TIME]);
        printf("        qnmpsi                 %12.4f\n", timings[QNMPSI_TIME]);
        printf("        theta * phi            %12.4f\n", timings[NONRES_TIME]);
        printf("        pulay mix phi          %12.4f\n", timings[MIXPSI_TIME]);
        printf("        H|phi>                 %12.4f\n", timings[HPSI_TIME]);


        printf("\n\n");
        printf("      mat & diag               %12.4f\n", timings[MATDIAG_TIME]);
        printf("-------------------------------------------\n");
        printf("        get_Hij                %12.4f\n", timings[HIJ_TIME]);
        printf("        get_matB               %12.4f\n", timings[MATB_TIME]);
        printf("        diagonalization        %12.4f\n", timings[DIAG_TIME]);
        printf("           orbit_dot_orbit     %12.4f\n", timings[ORBIT_DOT_ORBIT]);
        printf("           dot_product         %12.4f\n", timings[DOT_PRODUCT]);
        printf("           get_matB_qnm        %12.4f\n", timings[matB_QNM]);



        printf("\n\n");
        printf("      get_new_rho              %12.4f\n", timings[GET_NEW_RHO]);
        printf("-------------------------------------------\n");
        printf("        rho phi*phi            %12.4f\n", timings[RHO_PHI_TIME]);
        printf("        rho global sum         %12.4f\n", timings[RHO_SUM_TIME]);
        printf("        rho coarse to fine     %12.4f\n", timings[RHO_CTOF_TIME]);
        printf("        rho aug. term          %12.4f\n", timings[RHO_AUG_TIME]);
        printf("             rho_Qnm_mat       %12.4f\n", timings[RHO_QNM_MAT]);


        printf("\n\n");
        printf("    force                      %12.4f\n", timings[FORCE_TIME]);
        printf("-------------------------------------------\n");
        printf("      iiforce                  %12.4f\n", timings[IIFORCE_TIME]);
        printf("      lforce                   %12.4f\n", timings[LFORCE_TIME]);
        printf("      nlforce                  %12.4f\n", timings[NLFORCE_TIME]);
        printf("      nlccforce                %12.4f\n", timings[NLCCFORCE_TIME]);


        printf("\n\n");
        printf("      nlforce                  %12.4f\n", timings[NLFORCE_TIME]);
        printf("-------------------------------------------\n");
        printf("        get_all_partial_kbpsi  %12.4f\n", timings[PARTIAL_KBPSI]);
        printf("        rho_nm_mat             %12.4f\n", timings[RHO_NM_MAT]);
        printf("        partial_mat_nm_R       %12.4f\n", timings[PARTIAL_MAT_NM]);
        printf("        nlforce_par_Q          %12.4f\n", timings[NLFORCE_PAR_Q]);
        printf("        nlforce_par_rho        %12.4f\n", timings[NLFORCE_PAR_RHO]);
        printf("        nlforce_par_omega      %12.4f\n", timings[NLFORCE_PAR_OMEGA]);
        fflush(NULL);
    }                           /* end if */
    my_barrier();
    if (pct.gridpe == 0)
        printf("\n run done... \n");
}                               /* end write_timings */

/******/
