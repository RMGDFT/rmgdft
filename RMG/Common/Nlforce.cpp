/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "prototypes.h"



/*Set this to 1 to write out true NL force and the part
 * that comes from eigenvalues*/
#define VERBOSE 1

template void Nlforce<double> (double *, Kpoint<double> **Kptr);
template void Nlforce<std::complex<double> > (double * , Kpoint<std::complex<double>> **Kptr);

template <typename OrbitalType> void Nlforce (double * veff, Kpoint<OrbitalType> **Kptr)
{
    int ion, isp, index, gion, nion;
    int nh, size, size1;
    double *gamma, *par_gamma, *par_omega;
    SPECIES *sp;
    ION *iptr;
    int num_ions;
    fftw_plan p2;
    std::complex<double> *in, *out;
    double *newsintR_x, *newsintR_y, *newsintR_z, *qforce;
    double *newsintI_x, *newsintI_y, *newsintI_z, *tmp_force_gamma, *tmp_force_omega;
    int fpt0;
#if VERBOSE
    double *old_force, sum1x, sum1y, sum1z, sum2x, sum2y, sum2z;
    old_force = new double[ 3 * ct.num_ions];
#endif


    size1 = ct.num_kpts * ct.num_states * ct.num_ions * ct.max_nl;
    fpt0 = ct.fpt[0];


    num_ions = ct.num_ions;

    newsintR_x = new double[3 * size1];
    newsintR_y = newsintR_x + size1;
    newsintR_z = newsintR_y + size1;

    newsintI_x = new double[3 * size1];
    newsintI_y = newsintI_x + size1;
    newsintI_z = newsintI_y + size1;


    /*Initialize */
    for (isp = 0; isp < 3 * size1; isp++)
        newsintR_x[isp] = 0.0;
    for (isp = 0; isp < 3 * size1; isp++)
        newsintI_x[isp] = 0.0;


    /*Array for q-force */
    qforce = new double[ 3 * num_ions];
    tmp_force_gamma = new double[ 3 * num_ions];
    tmp_force_omega = new double[ 3 * num_ions];

    for (isp = 0; isp < 3 * num_ions; isp++)
    {
        qforce[isp] = 0.0;
        tmp_force_gamma[isp] = 0.0;
        tmp_force_omega[isp] = 0.0;
    }


    /*max for nh * (nh + 1) / 2 */
    size = (ct.max_nl + 1) * ct.max_nl / 2;
    
    gamma = new double[ size];
    par_gamma = new double[ 6 * size];
    par_omega = par_gamma + 3 * size;



    /*Loop over ions to setup newsintR_*, this loop is done separately to 
     * insure proper paralelization*/
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {
        /*Actual index of the ion under consideration*/
        gion = pct.nonloc_ions_list[ion];
        
        iptr = &ct.ions[gion];


        if (pct.idxptrlen[gion])
        {
            sp = &ct.sp[iptr->species];

#if !FDIFF_BETA
            in = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
            out = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
            p2 = fftw_plan_dft_3d (sp->nldim, sp->nldim, sp->nldim, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out),
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
#else
            p2 = NULL;
#endif


            PartialBetaxpsi (gion, p2, newsintR_x, newsintR_y, newsintR_z, newsintI_x, newsintI_y,
                              newsintI_z, iptr, Kptr);

            /*Release memery for plans */
#if !FDIFF_BETA
            fftw_destroy_plan (p2);
            fftw_free(out);
            fftw_free(in);
#endif




        }                       /*end if (pct.idxptrlen[ion]) */

        nh = ct.sp[iptr->species].nh;

        GetGamma (gamma, ion, nh, Kptr);
        nlforce_par_Q (veff, gamma, gion, iptr, nh, &qforce[3 * gion]);

    }                           /*end for(ion=0; ion<ions_max; ion++) */



    global_sums (newsintR_x, &size1, pct.grid_comm);
    global_sums (newsintR_y, &size1, pct.grid_comm);
    global_sums (newsintR_z, &size1, pct.grid_comm);

    global_sums (newsintI_x, &size1, pct.grid_comm);
    global_sums (newsintI_y, &size1, pct.grid_comm);
    global_sums (newsintI_z, &size1, pct.grid_comm);

    size1 = 3 * num_ions;
    global_sums (qforce, &size1, pct.img_comm);
        
    
    /*Add force calculated in nlforce1_par_Q */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];

        index = 3 * (ion);

#if VERBOSE
        old_force[index] = iptr->force[fpt0][0];
        old_force[index + 1] = iptr->force[fpt0][1];
        old_force[index + 2] = iptr->force[fpt0][2];
#endif
        iptr->force[fpt0][0] += get_vel_f() * qforce[index];
        iptr->force[fpt0][1] += get_vel_f() * qforce[index + 1];
        iptr->force[fpt0][2] += get_vel_f() * qforce[index + 2];
    }

#if VERBOSE
    printf("\n  Qforce");
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];

        index = 3 * (ion);
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f",
                    ion,
                    iptr->force[fpt0][0] - old_force[3 * ion],
                    iptr->force[fpt0][1] - old_force[3 * ion + 1],
                    iptr->force[fpt0][2] - old_force[3 * ion + 2]);
    }

    sum1x = 0.0;
    sum1y = 0.0;
    sum1z = 0.0;
    sum2x = 0.0;
    sum2y = 0.0;
    sum2z = 0.0;
#endif

    /*Loop over ions again */
    nion = -1;
    for (ion = 0; ion < pct.num_owned_ions; ion++)
    {
        /*Global index of owned ion*/
	gion = pct.owned_ions_list[ion];
	
        /* Figure out index of owned ion in nonloc_ions_list array, store it in nion*/
	do {
	    
	    nion++;
	    if (nion >= pct.num_nonloc_ions)
        {
            printf("\n Could not find matching entry in pct.nonloc_ions_list for owned ion %d", gion);
		rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in pct.nonloc_ions_list for owned ion ");
        }
	
	} while (pct.nonloc_ions_list[nion] != gion);
        
        iptr = &ct.ions[gion];


        nh = ct.sp[iptr->species].nh;

        /*partial_gamma(ion,par_gamma,par_omega, iptr, nh, p1, p2); */
        PartialGamma (gion, par_gamma, par_omega, nion, nh, newsintR_x, newsintR_y, newsintR_z,
                       newsintI_x, newsintI_y, newsintI_z, Kptr);
        nlforce_par_gamma (par_gamma, gion, nh, &tmp_force_gamma[3*gion]);


        nlforce_par_omega (par_omega, gion, nh, &tmp_force_omega[3*gion]);

    }                           /*end for(ion=0; ion<num_ions; ion++) */
    
    size1 = 3 * num_ions;
    global_sums (tmp_force_gamma, &size1, pct.img_comm);
    global_sums (tmp_force_omega, &size1, pct.img_comm);
    
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];

        index = 3 * (ion);
        iptr->force[fpt0][0] += tmp_force_gamma[index];
        iptr->force[fpt0][1] += tmp_force_gamma[index + 1];
        iptr->force[fpt0][2] += tmp_force_gamma[index + 2];
    }

    

#if VERBOSE
    if (pct.imgpe == 0)
    {
        printf ("\n\n True Non-local forces:");

        for (ion = 0; ion < ct.num_ions; ion++)
        {

            iptr = &ct.ions[ion];
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f",
                    ion,
                    iptr->force[fpt0][0] - old_force[3 * ion],
                    iptr->force[fpt0][1] - old_force[3 * ion + 1],
                    iptr->force[fpt0][2] - old_force[3 * ion + 2]);

            sum1x += iptr->force[fpt0][0] - old_force[3 * ion];
            sum1y += iptr->force[fpt0][1] - old_force[3 * ion + 1];
            sum1z += iptr->force[fpt0][2] - old_force[3 * ion + 2];
        }
        printf ("\n True NL sum in x, y and z directions: %e %e %e", sum1x, sum1y, sum1z);
    }

#endif


    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];

        index = 3 * (ion);

#if VERBOSE
        old_force[index] = iptr->force[fpt0][0];
        old_force[index + 1] = iptr->force[fpt0][1];
        old_force[index + 2] = iptr->force[fpt0][2];
#endif

        iptr->force[fpt0][0] += tmp_force_omega[index];
        iptr->force[fpt0][1] += tmp_force_omega[index + 1];
        iptr->force[fpt0][2] += tmp_force_omega[index + 2];
    }


#if VERBOSE
    if (pct.imgpe == 0)
    {
        printf ("\n\n Eigenvalue force:");

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
        printf ("\n Eigenvalue force sum in x, y and z directions: %e %e %e", sum2x, sum2y, sum2z);
    }
#endif

    delete[] par_gamma;
    delete[] gamma;
    delete[] tmp_force_gamma;
    delete[] tmp_force_omega;
    delete[] qforce;
    delete[] newsintR_x;
    delete[] newsintI_x;

#if VERBOSE
    delete[] old_force;
#endif

}
