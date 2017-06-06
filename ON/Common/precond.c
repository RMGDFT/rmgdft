/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include	"main.h"
#include	"prototypes_on.h"
#include "init_var.h"


void precond(double *x)
{
    double *work1, *work2;
    double *sp;
    double gamma;
    int size;
    int idx;
    int istate;

    int nx = ct.max_orbit_nx;

    my_malloc_init( work1, 4 * ct.max_orbit_size, double );
    my_malloc_init( work2, 4 * ct.max_orbit_size, double );

    gamma = get_gamma_precond(vtot_c, states[0].eig[0]);



    size = 0;
    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        sp = x + size;
        size += states[istate].size;
        int ixx = states[istate].ixmax - states[istate].ixmin + 1;
        int iyy = states[istate].iymax - states[istate].iymin + 1;
        int izz = states[istate].izmax - states[istate].izmin + 1;

        for (idx = 0; idx < states[istate].size; idx++)
        {
            work1[idx] = gamma * sp[idx];
        }
        /* compute the preconditioned steepest descent direction
         * -> work1 */

//        precond_mg_c(sp, work1, work2, istate);
        precond_mg(sp, work1, work2, istate);
          app_mask(istate, work1, 0);
//          ZeroBoundaryC(work1, ixx, iyy, izz);


        for (idx = 0; idx < states[istate].size; idx++)
        {
            sp[idx] = 0.5 * work1[idx];
        }

    }

    my_free(work1);
    my_free(work2);

}
