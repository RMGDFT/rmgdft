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

#ifndef CONST_H_INCLUDED
#define CONST_H_INCLUDED

/* AIX provides these (improperly) by default */
#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

/* Symbolic numeric constants */
#define		ZERO        0.0
#define		ONE         1.0
#define		TWO         2.0
#define		THREE       3.0
#define		FOUR        4.0
#define		FIVE        5.0
#define		SIX         6.0
#define		SEVEN       7.0
#define		EIGHT       8.0
#define		NINE        9.0
#define		TEN        10.0


#define		PI          3.14159265359
#define		twoPI       6.28318530718
#define		threePI     9.42477796077
#define		fourPI     12.56637061436


/* Physical constants and conversion factors obtained from NIST
 *
 *  http://physics.nist.gov/cgi-bin/cuu/Value?angstar|search_for=all!
 *
 */

/* electron charge in Coulomb */
#define     e_C         1.60217646e-19

/* Planck constant in SI */
#define     h_SI        6.626068e-34

/* unit conversions */
/* the numbers with higher precision were obtained by '1/const' */

/* conversions between hartree and eV */
#define     Ha_eV      27.2113845
#define     eV_Ha       0.0367493245336341

/* conversions between angstrom and bohr */
#define     a0_A        0.5291772108
#define     A_a0         1.88972612499359


/* conversion factor for atomic masses */
#define     mu_me    1822.88847984055
#define     me_mu       0.000548579910981427






/* Angular momentum states */
#define		S_STATE     0
#define		P_STATE     1
#define		D_STATE     2
#define		F_STATE     3
#define		G_STATE     4
#define		PX_STATE    1
#define		PY_STATE    2
#define		PZ_STATE    3


/* Neighbor list indices */
#define NB_N 0
#define NB_S 1
#define NB_E 2
#define NB_W 3
#define NB_U 4
#define NB_D 5
#define NB_SELF 7

/* relax method flags */
#define		FASTRELAX  0
#define		FIRE       1
#define		QUICK_MIN  2
#define         MD_MIN     3
#define		LBFGS      4

/* Molecular dynamics flags */
#define		MD_QUENCH       0
#define		MD_FASTRLX	1
#define		MD_CVE          2
#define		MD_CVT          3
#define		MD_CPT          4       /* not working yet */
#define		PLOT     	5
#define		PSIPLOT     	6
#define		BAND_STRUCTURE 	7
#define		NEB_RELAX       8

/* MD integration flags */
#define		ORDER_2         0
#define		ORDER_3         1
#define		ORDER_5         2

/* Temperature control flags */
#define		T_NOSE_CHAIN    0
#define		T_AND_SCALE     1

/* Boundary condition flags */
#define		PERIODIC        0
#define		CLUSTER         1
#define		SURFACE         2

/* Exchange-correlation flags */
#define       	LDA_PZ81        0
#define       	GGA_BLYP        1
#define       	GGA_XB_CP       2
#define       	GGA_XP_CP       3
#define       	GGA_PBE         4



/* Some stuff for timing and performance measurements */
#define TOTAL_TIME (0)
#define ORTHO_TIME (1)
//#define NL_TIME (2)
//#define NS_TIME (3)
#define EIG_TIME (4)
#define IMAGE_TIME (5)
//#define APPCIL_TIME (6)
//#define APPCIR_TIME (7)
#define RESTRICT_TIME (8)
#define EXPAND_TIME (9)
#define PACK_TIME (10)
#define INIT_TIME (11)
#define HARTREE_TIME (12)
#define DIAG_TIME (13)
#define LFORCE_TIME (14)
#define NLFORCE_TIME (15)
#define APPGRAD_TIME (16)
#define GATHER_TIME (17)
#define MG_EIGTIME (18)
#define INTERPOLATION_TIME (19)
#define RHO_TIME (20)
#define FORCE_TIME (21)
#define SCF_TIME (22)
//#define MD_LOOP (23)
#define NLCCFORCE_TIME (24)
#define IIFORCE_TIME (25)
#define MG_EIG_NLS_TIME (26)
#define MG_EIG_APPCIL_TIME (27)
#define MG_EIG_APPCIR_TIME (28)
#define MG_EIG_TRADE_TIME (29)
#define DIAG_NL_TIME (30)
#define DIAG_APPCIL_TIME (31)
#define DIAG_APPCIR_TIME (32)
#define GET_TE_XC_TIME (33)
#define GET_TE_II_TIME (34)
#define GET_TE_TIME (35)
#define SCF_XC_TIME (36)
#define INTERP_SETUP_TIME (37)
#define INTERP_EVAL_TIME (38)
#define DIAG_SUBDIAG1_TIME (39)
#define DIAG_MATRIX_TIME (40)
#define DIAG_WAVEUP_TIME (41)
#define DIAG_SUBDIAG1_LOOP_TIME (42)
#define DIAG_APP_A (43)
#define DIAG_APP_S (44)
#define DIAG_APP_B (45)
#define DIAG_DGEMM (46)
#define DIAG_GENVPSI_TIME (47)
#define DIAG_GLOB_SUMS (48)
#define SCF_BETAXPSI (49)
#define ALLOC_TIME (50)
#define ORTHO_NORM_PSI (52)
#define ORTHO_NEW_PSI (53)
#define ORTHO_GET_COEFF (54)
#define ORTHO_GLOB_SUM (55)
#define ORTHO_UPDATE_WAVES (56)
#define DIAG_APPCIR_TIME2 (57)
#define MG_EIG_GENVPSI_TIME (58)
#define MG_EIG_EIGVALUE_TIME (59)
#define MG_EIG_APPSMOOTH_TIME (60)
#define MG_EIG_MGRIDSOLV_TIME (61)
#define MG_EIG_PACK_TIME (62)
#define DIAG_NLS_TIME (63)
#define PREINIT_TIME (64)
#define FINISH_TIME (65)
#define DIAG_SCALAPACK_INIT (66)
#define DIAG_DISTMAT (67)
#define REAL_SUM_ALL_TIME (68)
#define GLOBAL_SUMS_TIME (69)
#define DIAG_BCAST_EIGS (70)
#define READ_PSEUDO_TIME 71
#define READ_CONTROL_TIME 72
#define ORTHO_GET_OVERLAPS 73
#define TRADE_MPI_TIME 74
#define ORTHO_CHOLESKY 75

#define LAST_TIME (80)


/* Occupation flags */
#define OCC_NONE 0
#define OCC_FD 1
#define OCC_GS 2
#define OCC_EF 3


/* Crystal lattice types */
/** Simple cubic lattice type.
 *  @doc Set input file value = 1 */
#define CUBIC_PRIMITIVE 	1

/** Face centered cubic lattice type. 
 *  @doc Set input file value = 2 */
#define CUBIC_FC		2

/** Bodycentered cubic lattice type. 
 *  @doc Set input file value = 3 */
#define CUBIC_BC		3

/** Hexagonal lattice type. 
 *  @doc Set input file value = 4 */
#define HEXAGONAL		4

#define TRIGONAL_PRIMITIVE	5
#define TETRAGONAL_PRIMITIVE	6
#define TETRAGONAL_BC           7

/** Orthorhombic lattice type. 
 *  @doc Set input file value = 8 */
#define ORTHORHOMBIC_PRIMITIVE  8

#define ORTHORHOMBIC_BASE_CENTRED 9
#define ORTHORHOMBIC_BC         10
#define ORTHORHOMBIC_FC 11
#define MONOCLINIC_PRIMITIVE 12
#define MONOCLINIC_BASE_CENTRED 13
#define TRICLINIC_PRIMITIVE 14

/* The real or imaginary part of a wavefunction */
#define PSI_REAL     0
#define PSI_IMAG     1

/* Type of async request passed to the mpi trade_images manager */
#define RMG_MPI_ISEND 1
#define RMG_MPI_IRECV 2


/* Types of finite differencing */
#define FULL_FD 1
#define CENTRAL_FD 2

/* Order of finite differencing for driver routines */
#define APP_CI_FOURTH 4
#define APP_CI_SIXTH 6

/* subspace diagonalization types */
#define SUBDIAG_LAPACK 0
#define SUBDIAG_SCALAPACK 1
#define SUBDIAG_MAGMA 2
#define SUBDIAG_ELPA 3

#endif /* CONST_H_INCLUDED */
