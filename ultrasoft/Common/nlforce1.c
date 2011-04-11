/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

/*Set this to 1 to write out true NL force and the part
 * that comes from eigenvalues*/
#define VERBOSE 1

void nlforce1 (REAL * veff)
{
    int ion, isp, index;
    int nh, size, size1;
    REAL *gamma, *par_gamma, *par_omega;
    SPECIES *sp;
    ION *iptr;
    int num_ions;
    fftwnd_plan p2;
    REAL *newsintR_x, *newsintR_y, *newsintR_z, *qforce;
    REAL *newsintI_x, *newsintI_y, *newsintI_z;
    int fpt0;
#if VERBOSE
    REAL *old_force, sum1x, sum1y, sum1z, sum2x, sum2y, sum2z;
    my_malloc (old_force, 3 * ct.num_ions, REAL);
#endif


    size1 = ct.num_kpts * ct.num_states * ct.num_ions * ct.max_nl;
    fpt0 = ct.fpt[0];


    num_ions = ct.num_ions;

    my_malloc (newsintR_x, 3 * size1, REAL);
    newsintR_y = newsintR_x + size1;
    newsintR_z = newsintR_y + size1;

#if !GAMMA_PT
    my_malloc (newsintI_x, 3 * size1, REAL);
    newsintI_y = newsintI_x + size1;
    newsintI_z = newsintI_y + size1;
#else
    newsintI_x = NULL;
    newsintI_y = NULL;
    newsintI_z = NULL;
#endif

    /*Initialize */
    for (isp = 0; isp < 3 * size1; isp++)
        newsintR_x[isp] = 0.0;
#if !GAMMA_PT
    for (isp = 0; isp < 3 * size1; isp++)
        newsintI_x[isp] = 0.0;
#endif


    /*Array for q-force */
    my_malloc (qforce, 3 * num_ions, REAL);

    for (isp = 0; isp < 3 * num_ions; isp++)
        qforce[isp] = 0.0;





    /*Loop over ions to setup newsintR_*, this loop is done separately to 
     * insure proper paralelization*/
    for (ion = 0; ion < num_ions; ion++)
    {
        iptr = &ct.ions[ion];


        if (pct.idxptrlen[ion])
        {
            sp = &ct.sp[iptr->species];

#if !FDIFF_BETA
            fftw_import_wisdom_from_string (sp->backward_wisdom);
            p2 = fftw3d_create_plan (sp->nldim, sp->nldim, sp->nldim,
                                     FFTW_BACKWARD, FFTW_USE_WISDOM);
            fftw_forget_wisdom ();
#else
            p2 = NULL;
#endif


            partial_betaxpsi (ion, p2, newsintR_x, newsintR_y, newsintR_z, newsintI_x, newsintI_y,
                              newsintI_z, iptr);

            /*Release memery for plans */
#if !FDIFF_BETA
            fftwnd_destroy_plan (p2);
#endif




        }                       /*end if (pct.idxptrlen[ion]) */

        nh = pct.prj_per_ion[ion];
        size = nh * (nh + 1) / 2;

        my_malloc (gamma, size, REAL);

        get_gamma (gamma, iptr, nh);
        nlforce1_par_Q (veff, gamma, ion, iptr, nh, &qforce[3 * ion]);

        my_free (gamma);


    }                           /*end for(ion=0; ion<ions_max; ion++) */



    global_sums (newsintR_x, &size1, pct.grid_comm);
    global_sums (newsintR_y, &size1, pct.grid_comm);
    global_sums (newsintR_z, &size1, pct.grid_comm);

#if !GAMMA_PT
    global_sums (newsintI_x, &size1, pct.grid_comm);
    global_sums (newsintI_y, &size1, pct.grid_comm);
    global_sums (newsintI_z, &size1, pct.grid_comm);
#endif

    size1 = 3 * num_ions;
    global_sums (qforce, &size1, pct.img_comm);

#if VERBOSE
    sum1x = 0.0;
    sum1y = 0.0;
    sum1z = 0.0;
    sum2x = 0.0;
    sum2y = 0.0;
    sum2z = 0.0;

    if (pct.imgpe == 0)
        printf ("\n\n True Non-local forces");
#endif

    /*Loop over ions again */
    for (ion = 0; ion < num_ions; ion++)
    {
        iptr = &ct.ions[ion];

#if VERBOSE
        old_force[3 * ion] = iptr->force[fpt0][0];
        old_force[3 * ion + 1] = iptr->force[fpt0][1];
        old_force[3 * ion + 2] = iptr->force[fpt0][2];
#endif


        nh = pct.prj_per_ion[ion];

        size = nh * (nh + 1) / 2;
        my_malloc (par_gamma, 6 * size, REAL);
        par_omega = par_gamma + 3 * size;

        /*partial_gamma(ion,par_gamma,par_omega, iptr, nh, p1, p2); */
        partial_gamma (ion, par_gamma, par_omega, iptr, nh, newsintR_x, newsintR_y, newsintR_z,
                       newsintI_x, newsintI_y, newsintI_z);
        nlforce1_par_gamma (par_gamma, ion, nh);

        /*Add force calculated in nlforce1_par_Q */
        index = 3 * (ion);
        iptr->force[fpt0][0] += ct.vel_f * qforce[index];
        iptr->force[fpt0][1] += ct.vel_f * qforce[index + 1];
        iptr->force[fpt0][2] += ct.vel_f * qforce[index + 2];

#if VERBOSE
        /*Print out true NL force */
        if (pct.imgpe == 0)
        {
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f", ion,
                    iptr->force[fpt0][0] - old_force[3 * ion],
                    iptr->force[fpt0][1] - old_force[3 * ion + 1],
                    iptr->force[fpt0][2] - old_force[3 * ion + 2]);

            sum1x += iptr->force[fpt0][0] - old_force[3 * ion];
            sum1y += iptr->force[fpt0][1] - old_force[3 * ion + 1];
            sum1z += iptr->force[fpt0][2] - old_force[3 * ion + 2];

            old_force[3 * ion] = iptr->force[fpt0][0];
            old_force[3 * ion + 1] = iptr->force[fpt0][1];
            old_force[3 * ion + 2] = iptr->force[fpt0][2];
        }
#endif


        nlforce1_par_omega (par_omega, ion, iptr, nh);

        my_free (par_gamma);
    }                           /*end for(ion=0; ion<num_ions; ion++) */


#if VERBOSE
    if (pct.imgpe == 0)
    {
        printf ("\n\n Eiegenvalue force:");

        for (ion = 0; ion < ct.num_ions; ion++)
        {

            iptr = &ct.ions[ion];
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f",
                    ion,
                    iptr->force[fpt0][0] - old_force[3 * ion],
                    iptr->force[fpt0][1] - old_force[3 * ion + 1],
                    iptr->force[fpt0][2] - old_force[3 * ion + 2]);

            sum2x += iptr->force[fpt0][0] - old_force[3 * ion];
            sum2y += iptr->force[fpt0][1] - old_force[3 * ion + 1];
            sum2z += iptr->force[fpt0][2] - old_force[3 * ion + 2];
        }
    }
#endif

    my_free (qforce);
    my_free (newsintR_x);
#if !GAMMA_PT
    my_free (newsintI_x);
#endif

#if VERBOSE

    if (pct.imgpe == 0)
        printf ("\n True NL sum in x, y and z directions: %e %e %e", sum1x, sum1y, sum1z);
    if (pct.imgpe == 0)
        printf ("\n Eigenvalue force sum in x, y and z directions: %e %e %e", sum2x, sum2y, sum2z);
    my_free (old_force);
#endif

}
