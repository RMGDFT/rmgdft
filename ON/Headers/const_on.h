/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/const.h *****
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
 *   
 * INPUTS
 *
 * OUTPUT
 *  
 * PARENTS
 *
 * CHILDREN
 * 
 * SEE ALSO
 *  
 * SOURCE
 */






#define INIT_FIREBALL  2
#define INIT_GAUSSIAN  3




/* Some stuff for timing and performance measurements */
#define NL_TIME (2)
#define MG_TIME (3)
#define APPCIL_TIME (5)
#define APPCIR_TIME (6)
#define MASK_TIME (17)
#define MATB_TIME (18)
#define PRECOND_TIME (19)
#define VXC_TIME (21)
#define SCF_TIME (22)
#define SX_TIME (23)
#define RES_TIME (24)
#define ORTHONORM_TIME (25)
#define COND_S_TIME (26)
#define GET_ORBIT (27)
#define OVERLAP_TIME (28)
#define CHOLESKY_TIME (29)
#define PSSYEV_TIME (30)
#define PSSYGST_TIME (31)
#define PSTRTRS_TIME (32)
#define LINE_MINTIME   (34)
#define ORBIT_DOT_ORBIT  (35)
#define ORBIT_DOT_ORBIT_H  (36)
#define DOT_PRODUCT (37)
#define GET_NEW_RHO (38)


#define MATDIAG_TIME (40)
#define UPDATEPOT_TIME  (41)
#define POTFC_TIME (42)
#define THETA_TIME (43)
#define QNMPSI_TIME (44)
#define NONRES_TIME (45)
#define HPSI_TIME (46)
#define DNMPSI_TIME (47)
#define MIXPSI_TIME (48)
#define HIJ_TIME (49)
#define INVMATB_TIME (50)
#define QUENCH_TIME (51)
#define RHO_PHI_TIME (52)
#define RHO_SUM_TIME (53)
#define RHO_CTOF_TIME (54)
#define RHO_AUG_TIME (55)
#define RHO_QNM_MAT (95)
#define matB_QNM (96)

#define PARTIAL_KBPSI (59)
#define RHO_NM_MAT (60)
#define PARTIAL_MAT_NM (61)
#define NLFORCE_PAR_Q (62)
#define NLFORCE_PAR_RHO (63)
#define NLFORCE_PAR_OMEGA (64)

#define READ_PSEUDO_TIME 71
#define READ_CONTROL_TIME 72




/******/
