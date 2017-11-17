/************************** SVN Revision Information **************************
 **    $Id: precond.c 4109 2017-06-06 21:00:32Z ebriggs $    **
******************************************************************************/
 
#include	"main.h"
#include	"prototypes_on.h"
#include "init_var.h"
#include <boost/circular_buffer.hpp>


void Precond(double *x)
{

    static boost::circular_buffer<double> rms(3);
    double rms_hist[3];
    static double beta=0.5;
    if(ct.scf_steps >= 3)
    {
        rms_hist[0] = rms.back();
        rms.pop_back();
        rms_hist[1] = rms.back();
        rms.pop_back();
        rms_hist[2] = rms.back();
        rms.pop_back();
        rms.push_back(rms_hist[2]);
        rms.push_back(rms_hist[1]);
        rms.push_back(rms_hist[0]);

        double ratio0 = rms_hist[0] / ct.rms;
        double ratio1 = rms_hist[1] / rms_hist[0];
        double ratio2 = rms_hist[2] / rms_hist[1];
#if 0
        if((ratio0 > ratio1) && (ratio1 > ratio2))
        {
             beta *= 2.0;
             if(ratio0 < 1.0) beta /= 4.0;
        }
        else if((ratio0 > ratio1) && (ratio1 < ratio2))
        {
             beta *= 2.0;
             if(ratio0 < 1.0) beta /= 4.0;
        }
        else if((ratio0 < ratio1) && (ratio1 > ratio2))
        { 
             if(ratio0 < 1.0) beta /= 2.0;
        }
        else if((ratio0 < ratio1) && (ratio1 < ratio2))
        {
            beta *= 2.0; 
             if(ratio0 < 1.0) beta /= 4.0;
        }
#endif
        if(ratio0 > ratio1)
            beta *= 1.5;
        else
            beta /= 1.5;

        if(beta > 1.5) beta = 1.5;
        printf("RATIO = %f  %f  %f  %f\n",ratio0, ratio1, ratio2, beta);
        printf("HIST  = %f  %f  %f  %f\n",ct.rms,rms_hist[0],rms_hist[1],rms_hist[2]);
    }
    rms.push_back(ct.rms);

    double *work1 = new double[4 * ct.max_orbit_size]();

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

        PrecondMg(psiR, work1, &states[istate]);

        for (int idx = 0; idx < states[istate].size; idx++)
        {
            psiR[idx] = beta*work1[idx];
        }

    }

    delete [] work1;

}
