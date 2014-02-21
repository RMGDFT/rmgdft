/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include	"main_on.h"



void precond(double *x)
{
    double *work1, *work2;
    double *sp;
    double gamma;
    int size;
    int idx;
    int istate;


    my_malloc_init( work1, 4 * ct.max_orbit_size, rmg_double_t );
    my_malloc_init( work2, 4 * ct.max_orbit_size, rmg_double_t );

    gamma = get_gamma_precond(vtot_c, states[0].eig[0]);

    size = 0;
    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        sp = x + size;
        size += states[istate].size;

        for (idx = 0; idx < states[istate].size; idx++)
        {
            work1[idx] = gamma * sp[idx];
        }
        /* compute the preconditioned steepest descent direction
         * -> work1 */
        precond_mg(sp, work1, work2, istate);
        app_mask(istate, work1, 0);

        for (idx = 0; idx < states[istate].size; idx++)
        {
            sp[idx] = work1[idx];
        }

    }

    my_free(work1);
    my_free(work2);

}
