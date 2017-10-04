/************************** SVN Revision Information **************************
 **    $Id: precond.c 4109 2017-06-06 21:00:32Z ebriggs $    **
******************************************************************************/
 
#include	"main.h"
#include	"prototypes_on.h"
#include "init_var.h"


void Precond(double *x)
{

    double *work1 = new double[4 * ct.max_orbit_size];
    double *work2 = new double[4 * ct.max_orbit_size];

    double gamma = get_gamma_precond(vtot_c, states[0].eig[0]);
    int size = 0;
    for (int istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        double *psiR = x + size;
        size += states[istate].size;
        int ixx = states[istate].ixmax - states[istate].ixmin + 1;
        int iyy = states[istate].iymax - states[istate].iymin + 1;
        int izz = states[istate].izmax - states[istate].izmin + 1;

        for (int idx = 0; idx < states[istate].size; idx++)
        {
            work1[idx] = gamma * psiR[idx];
        }
        /* compute the preconditioned steepest descent direction
         * -> work1 */

        ZeroBoundary(psiR, ixx, iyy, izz);
        PrecondMg(psiR, work1, &states[istate]);
        for (int idx = 0; idx < states[istate].size; idx++)
        {
            psiR[idx] = 0.5 * work1[idx];
        }
        ZeroBoundary(psiR, ixx, iyy, izz);

    }

    delete [] work2;
    delete [] work1;

}
