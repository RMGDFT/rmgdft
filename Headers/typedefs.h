/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*** QMD-MGDFT/main.h *****
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

#ifndef TYPEDEFS_H
#define TYPEDEFS_H 1

#if GPU_ENABLED
#include <cuda.h>
#include <cublas_v2.h>
#include <cublasXt.h>
#endif

#include <stdbool.h>
#include "fftw3.h"
#include "mpi.h"
#include "my_scalapack.h"


#include "pe_control.h"

/**@name STATE
 *
 * @memo Wavefunction storage structure */
typedef struct
{

    /** First iteration flag */
    int firstflag;

    /** Current estimate of the eigenvalue for this orbital (state). */
    double eig[2];

    /** Previous estimate */
    double oldeig[2];

    /** Wavefunction residual error computed by multigrid solver */
    double res;

    /** Points to the storage area for the real part of the orbital */
    double *psiR;
    /** Points to the storage area for the imaginary part of the orbital */
    double *psiI;


    /** Nuclear potential */
    double *vnuc;
    /** Hartree potential */
    double *vh;
    /** Exchange correlation potential */
    double *vxc;
    /** Total potential */
    double *vtot;

    /** Core charge for non-linear core corrections */
    double *rhocore;

    /** Total basis size on each processor (dimx*dimy*dimz) */
    int pbasis;

    /* Total basis size in a smoothing grid on each processor (dimx+2)*(dimy+2)*(dimz+2) */
    int sbasis;



    /** Volume element associated with each real space grid point */
    double vel;


    /** Occupation of the orbital */
    double occupation[2];
//    double oldeig;

    /* The offsets and the sizes of the grid that the orbital
     * is defined on relative to the global grid. These will
     * be used in the future for cluster boundary condition or
     * localized orbitals in an Order(N) formulation.
     */
    int xoff, yoff, zoff;
    int xsize, ysize, zsize;

    /** Index showing which k-point this orbital is associated with */
    int kidx;




    /* The ion on which this state is localized */
    int inum;

    /* index for functions with same localization */
    int loc_index;

    /* Actual Physical coordinates at current time step */
    int pe;
    double crds[3];
    double radius;
    int movable;
    int frozen;
    int index;

    int ixmin;
    int ixmax;
    int iymin;
    int iymax;
    int izmin;
    int izmax;
    int xfold;
    int yfold;
    int zfold;
    int ixstart;
    int iystart;
    int izstart;
    int ixend;
    int iyend;
    int izend;
    int orbit_nx, orbit_ny, orbit_nz;
    int size;
    /* Localization mask */
    char *lmask[4];

    int atomic_orbital_index;


    int n_orbital_same_center;
    int gaussian_orbital_index;

    /*8 Index of the orbital */
    int istate;
    int whichblock;
    int istate_in_block;
    int local_index;


    int ixmin_old;
    int ixmax_old;
    int iymin_old;
    int iymax_old;

    int atom_index;


} STATE;


#include "species.h"

/* Structure for storing species information for internal pseudopotentials */
typedef struct
{
    char name[4];
    double valence;
    double mass;
    double rc;
    int nlccflag;
    int maxl;
    int local;
} ISPECIES;




/* Nose control structure */
typedef struct
{

    /* number of atoms allowed to move */
    int N;

    /* ionic target temperature in Kelvin */
    double temp;

    /* ionic target kinetic energy */
    double k0;

    /* randomize velocity flag */
    bool randomvel;

    /* Nose oscillation frequency */
    double fNose;

    /* number of thermostats used */
    int m;

    /* thermostat positions,velocities,masses and forces */
    double xx[10];
    double xv[10];
    double xq[10];
    double xf[4][10];

} FINITE_T_PARM;


/** @name KPOINT
 * @memo Holds data specific to individual k-points.
 */
typedef struct
{

    /** The index of the k-point for backreferencing */
    int kidx;

    /** The k-point */
    double kpt[3];

    /** The corresponding vector */
    double kvec[3];

    /** The weight associated with the k-point */
    double kweight;

    /** The magnitude of the k-vector */
    double kmag;

    /* The orbital structure for this k-point */
    STATE *kstate;


    /* Mean min, and max wavefunction residuals for occupied space */
    double meanres;
    double minres;
    double maxres;

    /* Total energies */
    double ES;
    double NUC;
    double KE;
    double XC;
    double NL;
    double II;
    double TOTAL;

} KPOINT;


#include "rmg_control.h"


struct ION_ORBIT_OVERLAP
{
    short int xlow1;
    short int xhigh1;
    short int xlow2;
    short int xhigh2;
    short int xshift;

    short int ylow1;
    short int yhigh1;
    short int ylow2;
    short int yhigh2;
    short int yshift;

    short int zlow1;
    short int zhigh1;
    short int zlow2;
    short int zhigh2;
    short int zshift;

    short int flag;
};
typedef struct ION_ORBIT_OVERLAP ION_ORBIT_OVERLAP;


struct ORBIT_ORBIT_OVERLAP
{
    short int xlow1;
    short int xhigh1;
    short int xlow2;
    short int xhigh2;
    short int xshift;

    short int ylow1;
    short int yhigh1;
    short int ylow2;
    short int yhigh2;
    short int yshift;

    short int zlow1;
    short int zhigh1;
    short int zlow2;
    short int zhigh2;
    short int zshift;

    short int flag;
};
typedef struct ORBIT_ORBIT_OVERLAP ORBIT_ORBIT_OVERLAP;


#endif
