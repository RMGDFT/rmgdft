/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/***** Common-MGDFT/constrain.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
*                      Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void constraint(ION *)
 *   Enforces per atom constraints as set in the ct.ions[ion].constraint 3-vector.
 *   Atomic forces are projected onto the constraint vector in one of three ways
 *   based on the ct.ions[ion].contraint_type integer.
 *   constraint_type == 0::do nothing,
 *   constraint_type == 1::restrict to plane as defined by the constraint vector as normal,
 *   constraint_type == 2::restrict to direction as defined by the constraint vector.
 * INPUTS
 *   pointer to the ion that needs to have its constraint enforced.
 * OUTPUT
 *   no output
 * PARENTS
 *   rmg_fastrlx.c
 * CHILDREN
 *   
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include "main.h"

void constrain (ION * iptr)
{
    int fpt = ct.fpt[0];

    REAL fdotc =
        iptr->force[fpt][0] * iptr->constraint[0] +
        iptr->force[fpt][1] * iptr->constraint[1] +
        iptr->force[fpt][2] * iptr->constraint[2];

    switch (iptr->constraint_type)
    {
    case 2:                    /* Along vector, 1-D restriction (non-standard) */
        iptr->force[fpt][0] = fdotc * iptr->constraint[0];
        iptr->force[fpt][1] = fdotc * iptr->constraint[1];
        iptr->force[fpt][2] = fdotc * iptr->constraint[2];
        break;
    case 1:                    /* In plane, 2-D restriction (typical). */
        iptr->force[fpt][0] -= fdotc * iptr->constraint[0];
        iptr->force[fpt][1] -= fdotc * iptr->constraint[1];
        iptr->force[fpt][2] -= fdotc * iptr->constraint[2];
        break;
    case 0:                    /* Disable constraint */
        break;
    }
}
