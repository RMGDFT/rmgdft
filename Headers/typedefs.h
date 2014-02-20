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
#endif

#include "fftw.h"
#include "mpi.h"
#include "my_scalapack.h"

#if MPI
typedef struct
{
    int count;
} QMD_sem_t;
#endif


#include "pe_control.h"

/**@name STATE
 *
 * @memo Wavefunction storage structure */
typedef struct
{

    /** First iteration flag */
    int firstflag;

    /** Current estimate of the eigenvalue for this orbital (state). */
    rmg_double_t eig[2];

    /** Previous estimate */
    rmg_double_t oldeig[2];

    /** Wavefunction residual error computed by multigrid solver */
    rmg_double_t res;

    /** Points to the storage area for the real part of the orbital */
    rmg_double_t *psiR;
    /** Points to the storage area for the imaginary part of the orbital */
    rmg_double_t *psiI;


    /** Nuclear potential */
    rmg_double_t *vnuc;
    /** Hartree potential */
    rmg_double_t *vh;
    /** Exchange correlation potential */
    rmg_double_t *vxc;
    /** Total potential */
    rmg_double_t *vtot;

    /** dvhxc */
    rmg_double_t *dvhxc;

    /** Core charge for non-linear core corrections */
    rmg_double_t *rhocore;

    /** Grid dimension in the x-coordinate direction on this processor */
    int dimx;
    /** Grid dimension in the y-coordinate direction on this processor */
    int dimy;
    /** Grid dimension in the z-coordinate direction on this processor */
    int dimz;


    /** Grid spacings */
    rmg_double_t hxgrid;
    rmg_double_t hygrid;
    rmg_double_t hzgrid;


    /** Total basis size on each processor (dimx*dimy*dimz) */
    int pbasis;

    /* Total basis size in a smoothing grid on each processor (dimx+2)*(dimy+2)*(dimz+2) */
    int sbasis;

    /*8 Index of the orbital */
    int istate;


    /** Volume element associated with each real space grid point */
    rmg_double_t vel;


    /** Occupation of the orbital */
    rmg_double_t occupation[2];
//    rmg_double_t oldeig;

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
    rmg_double_t crds[3];
    rmg_double_t radius;
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

} STATE;


#include "species.h"

/* Structure for storing species information for internal pseudopotentials */
typedef struct
{
    char name[4];
    rmg_double_t valence;
    rmg_double_t mass;
    rmg_double_t rc;
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
    rmg_double_t temp;

    /* ionic target kinetic energy */
    rmg_double_t k0;

    /* randomize velocity flag */
    int randomvel;

    /* Nose oscillation frequency */
    rmg_double_t fNose;

    /* number of thermostats used */
    int m;

    /* thermostat positions,velocities,masses and forces */
    rmg_double_t xx[10];
    rmg_double_t xv[10];
    rmg_double_t xq[10];
    rmg_double_t xf[4][10];

} FINITE_T_PARM;


/** @name KPOINT
 * @memo Holds data specific to individual k-points.
 */
typedef struct
{

    /** The index of the k-point for backreferencing */
    int kidx;

    /** The k-point */
    rmg_double_t kpt[3];

    /** The corresponding vector */
    rmg_double_t kvec[3];

    /** The weight associated with the k-point */
    rmg_double_t kweight;

    /** The magnitude of the k-vector */
    rmg_double_t kmag;

    /* The orbital structure for this k-point */
    STATE *kstate;


    /* Mean min, and max wavefunction residuals for occupied space */
    rmg_double_t meanres;
    rmg_double_t minres;
    rmg_double_t maxres;

    /* Total energies */
    rmg_double_t ES;
    rmg_double_t NUC;
    rmg_double_t KE;
    rmg_double_t XC;
    rmg_double_t NL;
    rmg_double_t II;
    rmg_double_t TOTAL;

} KPOINT;

#define         MAX_STATES  3000

#include "rmg_control.h"

/* Extern declaration for the processor control structure */
extern PE_CONTROL pct;

#endif
