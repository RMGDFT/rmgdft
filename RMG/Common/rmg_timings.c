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
 *   print the timings in the standard output.
 * PARENTS
 *   too many
 * CHILDREN
 *   no
 * SOURCE
 */


#include "main.h"
#include <time.h>
#include <stdio.h>

#if HYBRID_MODEL
#include "hybrid.h"
#endif

#if PAPI_PERFMON
#undef kill
#include <papi.h>
#endif




REAL timings[LAST_TIME];

// For the hybrid model case we need to lock the timings array with a mutex and do some
// adjustments to the timings to account for multiple threads so the function is defined
// in hybrid.c instead of here.
#if !HYBRID_MODEL
void rmg_timings (int what, REAL time)
{

    timings[what] += time;
}                               /* end rmg_timings */
#endif




#include <sys/time.h>

REAL my_crtc (void)
{
    struct timeval t1;
    gettimeofday (&t1, NULL);
    return t1.tv_sec + 1e-6 * t1.tv_usec;
}


#define printf_timing_line1( _title_, _index_ ) printf(_title_ "%12.2f (%5.1f%%)\n", \
                                                       timings[_index_],\
                                                       timings[_index_] / total_time *100.0 );


#define printf_timing_line2( _title_, _index_ ) printf(_title_ "%12.2f (%5.1f%%)    %12.2f\n", \
                                                       timings[_index_],\
                                                       timings[_index_] / total_time *100.0, \
                                                       timings[_index_] / md_steps)


#define printf_timing_line3( _title_, _index_ ) printf(_title_ "%12.2f (%5.1f%%)    %12.2f     %12.2f\n", \
                                                       timings[_index_],\
                                                       timings[_index_] / total_time *100.0, \
                                                       timings[_index_] / md_steps, \
                                                       timings[_index_] / total_scf_steps)



/* Outputs timing information */
void write_timings (void)
{
    REAL total_time, FLOPS, TOTAL_FLOPS=0.0;
    int i, md_steps, total_scf_steps, ithread;

#if PAPI_PERFMON
    long long Papi_values[4];
    Papi_values[0] = 0;
    Papi_values[1] = 0;
    Papi_values[2] = 0;
    Papi_values[3] = 0;

#if HYBRID_MODEL
    enter_threaded_region();
    scf_barrier_init(THREADS_PER_NODE);
    for(ithread = 0;ithread < THREADS_PER_NODE;ithread++) {
          thread_control[ithread].job = HYBRID_FINALIZE_PAPI;
    }

    // Thread tasks are set up so wake them
    wake_threads(THREADS_PER_NODE);

    // Then wait for them to finish this task
    wait_for_threads(THREADS_PER_NODE);
    scf_barrier_destroy();
    leave_threaded_region();

#pragma omp parallel for
    for(i = 0;i < THREADS_PER_NODE;i++)
        Papi_finalize_omp_threads(i);

    printf("\n\nPer thread performance data for 1 PE\n");
    for(ithread = 0;ithread < THREADS_PER_NODE;ithread++) {
        FLOPS = ((double)Papi_thread_flops(ithread)) / (my_crtc () - ct.time0);
        printf("  MFLOPS for Posix thread %d  = %f \n", ithread, FLOPS / 1.0e6);
        TOTAL_FLOPS += Papi_thread_flops(ithread);
    }

    for(ithread = 0;ithread < THREADS_PER_NODE;ithread++) {
        FLOPS = ((double)ct.OpenMpFlopCount[ithread]) / (my_crtc () - ct.time0);
        printf("  MFLOPS for OpenMP thread %d = %f \n", ithread, FLOPS / 1.0e6);
        TOTAL_FLOPS += ct.OpenMpFlopCount[ithread];
    }
#endif

    printf("Total performance for this PE = %f  MFLOPS\n\n", TOTAL_FLOPS / (my_crtc () - ct.time0) / 1.0e6);

#endif


    /*Since we count SCF and MD steps from 0, we are undercounting by one */
    md_steps = ct.md_steps + 1;
    total_scf_steps = ct.total_scf_steps + 1;

    /*Total time will be time from the start until now */
    total_time = my_crtc () - ct.time0;

    timings[TOTAL_TIME] = total_time;




    printf ("\n\n\n");
    printf ("------------------------------ TIMING INFORMATION ---------------------------\n");
    printf ("                          Whole Run             Per MD Step      Per SCF Step\n");

    printf_timing_line3 (" Total time            ", TOTAL_TIME);
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line1 ("   Pre-Initialization  ", PREINIT_TIME);
    printf_timing_line1 ("     Read Control File ", READ_CONTROL_TIME);
    printf_timing_line1 ("     Read Pseudo File  ", READ_PSEUDO_TIME);
    printf_timing_line1 ("   Pre-Initialization  ", PREINIT_TIME);
    printf_timing_line1 ("   Initialization      ", INIT_TIME);
    printf_timing_line3 ("   Scf steps           ", SCF_TIME);
    printf_timing_line3 ("   Get_te              ", GET_TE_TIME);
    printf_timing_line3 ("   Diagonalization     ", DIAG_TIME);
    printf_timing_line3 ("   Get/Release_mem     ", ALLOC_TIME);
    printf_timing_line2 ("   Force               ", FORCE_TIME);
    printf_timing_line1 ("   Finish up time      ", FINISH_TIME);

    printf ("\n");
    printf_timing_line3 (" Scf steps             ", SCF_TIME);
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line3 ("   Hartree             ", HARTREE_TIME);
    printf_timing_line3 ("   Orthogonalization   ", ORTHO_TIME);
    printf_timing_line3 ("   Get vxc             ", SCF_XC_TIME);
    printf_timing_line3 ("   Get rho             ", RHO_TIME);
    printf_timing_line3 ("   Mg_eig              ", EIG_TIME);
    printf_timing_line3 ("   Beta x psi          ", SCF_BETAXPSI);

#if MD_TIMERS
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line3 ("     App_nl & App_ns   ", MG_EIG_NLS_TIME);
    printf_timing_line3 ("     App_cil & smooth  ", MG_EIG_APPCIL_TIME);
    printf_timing_line3 ("     App_cir           ", MG_EIG_APPCIR_TIME);
    printf_timing_line3 ("     Trade_image       ", MG_EIG_TRADE_TIME);
    printf_timing_line3 ("     Genvpsi           ", MG_EIG_GENVPSI_TIME);
    printf_timing_line3 ("     Calculate Eigenval", MG_EIG_EIGVALUE_TIME);
    printf_timing_line3 ("     App_smooth        ", MG_EIG_APPSMOOTH_TIME);
    printf_timing_line3 ("     Mg_solv           ", MG_EIG_MGRIDSOLV_TIME);
    printf_timing_line3 ("     Packing           ", MG_EIG_PACK_TIME);
#endif


#if MD_TIMERS
    printf_timing_line3 ("   Orthogonalization   ", ORTHO_TIME);
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line3 ("     Normalization     ", ORTHO_NORM_PSI);
    printf_timing_line3 ("     Psi overlaps      ", ORTHO_GET_OVERLAPS);
    printf_timing_line3 ("     Calculate coeffs  ", ORTHO_GET_COEFF);
    printf_timing_line3 ("     Global sums       ", ORTHO_GLOB_SUM);
    printf_timing_line3 ("     Update psi        ", ORTHO_UPDATE_WAVES);
#endif


    printf ("\n");
    printf_timing_line3 (" Get rho               ", RHO_TIME);
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line3 ("   Interpolation       ", INTERPOLATION_TIME);

#if MD_TIMERS
    if (ct.interp_flag)
    {
        printf ("-----------------------------------------------------------------------------\n");
        printf_timing_line3 ("     B-Spline coeff    ", INTERP_SETUP_TIME);
        printf_timing_line3 ("     B-Spline evaluat  ", INTERP_EVAL_TIME);
    }
#endif



    printf ("\n");
    printf_timing_line3 (" Get_te                ", GET_TE_TIME);
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line3 ("   XC Energy           ", GET_TE_XC_TIME);
    printf_timing_line3 ("   Ion-Ion Energy      ", GET_TE_II_TIME);

    printf ("\n");
    printf_timing_line3 (" Diagonalization       ", DIAG_TIME);
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line3 ("   Global Sums         ", DIAG_GLOB_SUMS);
    printf_timing_line3 ("   Broadcast eigenvals ", DIAG_BCAST_EIGS);
#if GAMMA_PT
    printf_timing_line3 ("   ScaLapack init      ", DIAG_SCALAPACK_INIT);
    printf_timing_line3 ("   Matrix distribution ", DIAG_DISTMAT);
#endif
    printf_timing_line3 ("   ScaLapack operations", DIAG_MATRIX_TIME);
    printf_timing_line3 ("   Wavefunc update     ", DIAG_WAVEUP_TIME);
#if GAMMA_PT
    printf_timing_line3 ("   Apply A operator    ", DIAG_APP_A);
#if MD_TIMERS
    printf_timing_line3 ("    -App_nl            ", DIAG_NL_TIME);
    printf_timing_line3 ("    -Genvpsi           ", DIAG_GENVPSI_TIME);
    printf_timing_line3 ("    -App_cir           ", DIAG_APPCIR_TIME);
    printf_timing_line3 ("    -App_cil           ", DIAG_APPCIL_TIME);
#endif
    printf_timing_line3 ("   Apply S operator    ", DIAG_APP_S);
    printf_timing_line3 ("   Apply B operator    ", DIAG_APP_B);
#if MD_TIMERS
    printf_timing_line3 ("    -App_cir           ", DIAG_APPCIR_TIME2);
#endif
    printf_timing_line3 ("   Setup Aij, Bij, Cij ", DIAG_DGEMM);

    /*Non-gamma point */
#else
    printf_timing_line3 ("   Subdiag1            ", DIAG_SUBDIAG1_TIME);
#if MD_TIMERS
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line3 ("     App_nl & App_ns   ", DIAG_NLS_TIME);
    printf_timing_line3 ("     App_cil           ", DIAG_APPCIL_TIME);
    printf_timing_line3 ("     App_cir           ", DIAG_APPCIR_TIME);
    printf_timing_line3 ("     Setting A and B   ", DIAG_SUBDIAG1_LOOP_TIME);
#endif
#endif



    printf ("\n");
    printf_timing_line2 (" Force                 ", FORCE_TIME);
    printf ("-----------------------------------------------------------------------------\n");
    printf_timing_line2 ("   Non-local force     ", NLFORCE_TIME);
#if MD_TIMERS
    printf_timing_line2 ("   Local force         ", LFORCE_TIME);
    printf_timing_line2 ("   Non-linear core     ", NLCCFORCE_TIME);
    printf_timing_line2 ("   Ion-Ion force       ", IIFORCE_TIME);
#endif

#if MD_TIMERS
    printf ("\n");
    /*These are timers that count how long it took to do any trade_image (0,2,3,5) and 
     * gather and scatter*/
    printf_timing_line3 (" Trade images          ", IMAGE_TIME);
    printf_timing_line3 (" Trade images mpi time ", TRADE_MPI_TIME);
    printf_timing_line3 (" Real_sum_all          ", REAL_SUM_ALL_TIME);
    printf_timing_line3 (" GLOBAL_SUMS_TIME      ", GLOBAL_SUMS_TIME);
    printf_timing_line3 (" Gather                ", GATHER_TIME);
#endif

    if (ct.md_steps)
    {
        printf ("\n\n");
        printf ("%d SCF steps were done in %d MD steps\n", ct.total_scf_steps, ct.md_steps);
        printf (" Average: %.1f SCF steps per MD step\n\n",
                (REAL) ct.total_scf_steps / (REAL) ct.md_steps);
    }


}                               /* end write_timings */


/* Resets timers at the beginning of MD run, so that end timings are not "poluted"
 *  * by initial quench run*/
void reset_timers (void)
{
    int i;

    for (i = 0; i < LAST_TIME; i++)
    {
        if ((i == TOTAL_TIME) || (i == INIT_TIME))
            continue;
        timings[i] = 0.0;
    }

}                               /*end reset_timers */

/******/
