/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/***** RMG/Common/reinit_ionic_pp.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2010  : Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void reinit_ionic_pp (STATE * states, REAL * vnuc, REAL * rhocore, REAL * rhoc)
 *   drive routine for updating systemic potentials due to motion of ions.
 * INPUTS
 *   states: all wave functions (see main.h)
 *   vnuc: pseudopotential
 *   rhocore: charge density of core electrons, only useful when we 
 *            include non-linear core correction for pseudopotential.
 *   rhoc:    compensating charge density
 * OUTPUT
 *   all the above inputs are updated
 * PARENTS
 *   fastrlx.c, moldyn.c
 * CHILDREN
 *   init_nuc,get_QI,get_nlop,get_weight,get_qqq,betaxpsi,mix_betaxpsi
 * SOURCE
 */

#include "main.h"

static gpu_weight_alloc=0;

void reinit_ionic_pp (STATE * states, REAL * vnuc, REAL * rhocore, REAL * rhoc)
{

    /* Update items that change when the ionic coordinates change */
    init_nuc (vnuc, rhoc, rhocore);
    get_QI ();
    get_nlop ();

    /*Other things that need to be recalculated when ionic positions change */
    get_weight ();
    get_qqq ();

#if GPU_ENABLED
    // If gpu weight buffer has not been setup yet or size has changed must take care of allocation
    if((ct.gpu_weight == NULL) || (gpu_weight_alloc < (pct.P0_BASIS * pct.num_tot_proj * sizeof(REAL)))) {
        if(ct.gpu_weight != NULL) {
            cudaFree(ct.gpu_weight); 
            ct.gpu_weight = NULL;
        }
        if(pct.num_tot_proj) {
            if( cudaSuccess != cudaMalloc((void **)&ct.gpu_weight , pct.P0_BASIS * pct.num_tot_proj * sizeof(REAL) ))
                error_handler("cudaMalloc failed for: gpu_weight\n");
        }
        gpu_weight_alloc = pct.P0_BASIS * pct.num_tot_proj * sizeof(REAL);
    }
    // Transfer copy of weights to GPU
    cublasSetVector( pct.P0_BASIS * pct.num_tot_proj, sizeof( REAL ), pct.weight, 1, ct.gpu_weight, 1 );
#endif

    if((ct.gpu_Bweight == NULL) || (gpu_weight_alloc < (pct.P0_BASIS * pct.num_tot_proj * sizeof(REAL)))) {
        if(ct.gpu_Bweight != NULL) {
            cudaFree(ct.gpu_Bweight); 
            ct.gpu_Bweight = NULL;
        }
        if(pct.num_tot_proj) {
            if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Bweight , pct.P0_BASIS * pct.num_tot_proj * sizeof(REAL) ))
                error_handler("cudaMalloc failed for: gpu_Bweight\n");
        }
    }
    // Transfer copy of weights to GPU
    cublasSetVector( pct.P0_BASIS * pct.num_tot_proj, sizeof( REAL ), pct.Bweight, 1, ct.gpu_Bweight, 1 );
#endif

    if (!verify ("calculation_mode", "Band Structure Only"))
    {
        betaxpsi (states);
        mix_betaxpsi(0);
    }

}                               /* end reinit_ionic_pp */


/******/
