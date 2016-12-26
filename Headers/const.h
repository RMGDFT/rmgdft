/************************** SVN Revision Information **************************
 **    $Id: const.h 2072 2013-11-01 14:50:01Z ebriggs $    **
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

// Discretization types
#define         MEHRSTELLEN_DISCRETIZATION      0
#define         CENTRAL_DISCRETIZATION          1

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


#define		PI          3.14159265358979323
#define		twoPI       6.28318530717958646
#define		threePI     9.42477796076937969
#define		fourPI     12.56637061435917292
#define         SQRT2       1.41421356237309504880
#define         SQRT3       1.73205080756887729352


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
#define		TDDFT       10

/* MD integration flags */
#define		ORDER_2         0
#define		ORDER_3         1
#define		ORDER_5         2

/* Temperature control flags */
#define		T_NOSE_CHAIN    0
#define		T_AND_SCALE     1

/* Exchange-correlation flags */
#define       	LDA_PZ81        0
#define       	GGA_BLYP        1
#define       	GGA_XB_CP       2
#define       	GGA_XP_CP       3
#define       	GGA_PBE         4
#define       	MGGA_TB09       5
#define       	VDW             6
#define       	AUTO_XC         99


/******/

/* Occupation flags */
#define OCC_NONE 0
#define OCC_FD 1
#define OCC_GS 2
#define OCC_EF 3
#define OCC_MV 4
#define OCC_MP 5


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
#define PSI_double     0
#define PSI_IMAG     1

/* subspace diagonalization types */
#define SUBDIAG_LAPACK 0
#define SUBDIAG_SCALAPACK 1
#define SUBDIAG_MAGMA 2
#define SUBDIAG_AUTO 3
#define SUBDIAG_ELPA 6
#define SUBDIAG_MAGMAFS 20
#define SUBDIAG_LAPACKFS 21
#define SUBDIAG_SCALAPACKFS 22

// Run type flags
#define RANDOM_START 0
#define RESTART 1
#define LCAO_START 2
#define INIT_FIREBALL  3
#define INIT_GAUSSIAN  4
#define Start_TDDFT  5
#define Restart_TDDFT  6
#define MODIFIED_LCAO_START  7

// 3-d interpolation types
#define CUBIC_POLYNOMIAL_INTERPOLATION 0
#define BSPLINE_INTERPOLATION 1
#define PROLONG_INTERPOLATION 2
#define FFT_INTERPOLATION 3

// Interpolation type flags for pseudopotentials
#define LINEAR_INTERPOLATION 0
#define LOGARITHMIC_INTERPOLATION 1

// Kohn-sham solver types
#define MULTIGRID_SOLVER 0
#define DAVIDSON_SOLVER 1
#define POISSON_PFFT_SOLVER 1

// Fft filtering types
#define LOW_PASS 0
#define HIGH_PASS 1

// Force types
#define WAVEFUNCTION_DERIVATIVE 0
#define PROJECTOR_DERIVATIVE 1

//Charge Analysis Method
#define CHARGE_ANALYSIS_NONE 0 
#define CHARGE_ANALYSIS_VORONOI 1 

// Constants for InitLocalObject
#define ATOMIC_LOCAL_PP  0
#define ATOMIC_RHO       1
#define ATOMIC_RHOCOMP   2
#define ATOMIC_RHOCORE   3

#endif /* CONST_H_INCLUDED */
