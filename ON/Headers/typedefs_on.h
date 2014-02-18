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


/* multigrid-parameter structure */
typedef struct
{

    /* number of global-grid pre/post smoothings and timestep */
    double gl_step;
    int gl_pre;
    int gl_pst;

    /* timestep for the subiteration */
    double sb_step;

    /* timestep for the Richardson-Iteration */
    double ri_step;

    /* lowest MG level */
    int levels;

    /* Number of Mu-cycles to use */
    int mucycles;

    /* Number of Smoother iterations on the coarsest level */
    int coarsest_steps;

} MG_PARM;

