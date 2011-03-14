/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

   get_nlop.c


   Sets up the ket part of the non-local operators.



 */




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"
#include "fftw.h"

void init_derweight();
void init_weight();

void get_nlop(void)
{
    int ion, idx, ip, isp;
    int tot_prj, item, ion1, index;
    int PROJECTOR_SPACE, prjcount;
    double *beta, *beta_x, *beta_y, *beta_z;
    SPECIES *sp;
    ION *iptr;
    fftwnd_plan p1, p2;
    int overlap;
    int coarse_size, st1;
    double *fftw_phase_sin,*fftw_phase_cos;
    double vect[3], nlcrds[3];

    /*Pointer to the result of forward transform on the coarse grid */
    fftw_complex *fptr, *fptr_x, *fptr_y, *fptr_z;
    fftw_complex *beptr, *gbptr;



    double time2 = my_crtc();
    /*Do forward transform for each species and store results on the coarse grid */
    init_weight ();
    /*The same for derivative of beta */
    init_derweight ();

    if (pct.gridpe == 0)
        printf ("\n init: FFTW initialization finished, it took %.1f s", my_crtc () - time2);

    /*Get memory to store the phase factor applied to the forward Fourier transform
     *      * and to store the backwards transform*/
    my_malloc( beptr, 2 * ct.max_nlpoints, fftw_complex );
    if (beptr == NULL)
        error_handler ("can't allocate memory\n");

    gbptr = beptr + ct.max_nlpoints;

    
    my_malloc( fftw_phase_sin, 2 * ct.max_nlpoints, double );
    fftw_phase_cos = fftw_phase_sin + ct.max_nlpoints;



    /*
     * PROJECTOR_SPACE = ct.max_nlpoints * ct.max_nl;
     */

    my_barrier();

    /*  get total number of projectors on this processor */
    /*  pct.n_ion_center: number of ions whose nl projector overlap
     *  with the states on this processor */

    pct.n_ion_center = 0;
    tot_prj = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        overlap = 0;
        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        {
            index = (st1 - ct.state_begin) * ct.num_ions + ion;
            if (ion_orbit_overlap_region_nl[index].flag == 1)
                overlap = 1;
        }
        if (overlap == 1)
        {
            pct.ionidx[pct.n_ion_center] = ion;
            pct.n_ion_center += 1;
            tot_prj += ct.sp[ct.ions[ion].species].num_projectors;
        }
    }

    PROJECTOR_SPACE = ct.max_nlpoints * tot_prj;
    if (projectors != NULL)
        my_free(projectors);
    my_malloc_init( projectors, PROJECTOR_SPACE, REAL );

    /*allocate memorry for weight factor of partial_beta/partial_x */
    if (projectors_x != NULL)
        my_free(projectors_x);
    my_malloc_init( projectors_x, PROJECTOR_SPACE, REAL );

    /*allocate memorry for weight factor of partial_beta/partial_y */
    if (projectors_y != NULL)
        my_free(projectors_y);
    my_malloc_init( projectors_y, PROJECTOR_SPACE, REAL );

    /*allocate memorry for weight factor of partial_beta/partial_z */
    if (projectors_z != NULL)
        my_free(projectors_z);
    my_malloc_init( projectors_z, PROJECTOR_SPACE, REAL );


    for (isp = 0; isp < ct.num_species; isp++)
    {
        sp = &ct.sp[isp];

        p2 = fftw3d_create_plan(sp->nldim, sp->nldim, sp->nldim,
                FFTW_BACKWARD, FFTW_MEASURE | FFTW_USE_WISDOM);
        sp->backward_wisdom = fftw_export_wisdom_to_string();
        fftwnd_destroy_plan(p2);
        fftw_forget_wisdom();
    }

    for (ion = 0; ion < ct.num_ions; ion++)
        pct.prj_per_ion[ion] = ct.sp[ct.ions[ion].species].num_projectors;

    /* Loop over all the ions on this processor */

    beta = projectors;
    beta_x = projectors_x;
    beta_y = projectors_y;
    beta_z = projectors_z;



    prjcount = 0;
    for (ion1 = 0; ion1 < pct.n_ion_center; ion1++)
    {
        ion = pct.ionidx[ion1];
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        fftw_import_wisdom_from_string(sp->backward_wisdom);
        p2 = fftw3d_create_plan(sp->nldim, sp->nldim, sp->nldim, FFTW_BACKWARD, FFTW_USE_WISDOM);
        fftw_forget_wisdom();

        /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
        /*xtal vector between ion and left bottom corner of the box */

        vect[0] = iptr->xtal[0] - iptr->nlxcstart;
        vect[1] = iptr->xtal[1] - iptr->nlycstart;
        vect[2] = iptr->xtal[2] - iptr->nlzcstart;

        /*Substract vector between left bottom corner of the box and center of the box */
        vect[0] -= (sp->nldim / 2) / (REAL) ct.psi_nxgrid;
        vect[1] -= (sp->nldim / 2) / (REAL) ct.psi_nygrid;
        vect[2] -= (sp->nldim / 2) / (REAL) ct.psi_nzgrid;

        /*The vector we are looking for should be */
        to_cartesian (vect, nlcrds);


        coarse_size = sp->nldim *sp->nldim *sp->nldim ;
        
        /*Calculate the phase factor */
        find_phase (sp->nldim, nlcrds, fftw_phase_sin, fftw_phase_cos);

        /*Temporary pointer to the already calculated forward transform */
        fptr = sp->forward_beta;
        fptr_x = sp->forward_derbeta_x;
        fptr_y = sp->forward_derbeta_y;
        fptr_z = sp->forward_derbeta_z;

        /* Loop over radial projectors */
        for (ip = 0; ip < sp->num_projectors; ip++)
        {


            /****************** beta *************/
            /*Apply the phase factor   */
            for (idx = 0; idx < coarse_size; idx++)
            {
                gbptr[idx].re =
                    fptr[idx].re *fftw_phase_cos[idx] +
                    fptr[idx].im * fftw_phase_sin[idx];
                gbptr[idx].im =
                    fptr[idx].im * fftw_phase_cos[idx] -
                    fptr[idx].re * fftw_phase_sin[idx];
            }


            /*Do the backwards transform */
            fftwnd_one (p2, gbptr, beptr);
            /*This takes and stores the part of beta that is useful for this PE */
            assign_weight (sp, beptr, &beta[prjcount]);

            /****************** beta_X *************/
            /*Apply the phase factor   */
            for (idx = 0; idx < coarse_size; idx++)
            {
                gbptr[idx].re =
                    fptr_x[idx].re *fftw_phase_cos[idx] +
                    fptr_x[idx].im * fftw_phase_sin[idx];
                gbptr[idx].im =
                    fptr_x[idx].im * fftw_phase_cos[idx] -
                    fptr_x[idx].re * fftw_phase_sin[idx];
            }


            /*Do the backwards transform */
            fftwnd_one (p2, gbptr, beptr);
            /*This takes and stores the part of beta that is useful for this PE */
            assign_weight (sp, beptr, &beta_x[prjcount]);


            /****************** beta_Y *************/
            /*Apply the phase factor   */
            for (idx = 0; idx < coarse_size; idx++)
            {
                gbptr[idx].re =
                    fptr_y[idx].re *fftw_phase_cos[idx] +
                    fptr_y[idx].im * fftw_phase_sin[idx];
                gbptr[idx].im =
                    fptr_y[idx].im * fftw_phase_cos[idx] -
                    fptr_y[idx].re * fftw_phase_sin[idx];
            }


            /*Do the backwards transform */
            fftwnd_one (p2, gbptr, beptr);
            /*This takes and stores the part of beta that is useful for this PE */
            assign_weight (sp, beptr, &beta_y[prjcount]);



            /****************** beta_Z *************/
            /*Apply the phase factor   */
            for (idx = 0; idx < coarse_size; idx++)
            {
                gbptr[idx].re =
                    fptr_z[idx].re *fftw_phase_cos[idx] +
                    fptr_z[idx].im * fftw_phase_sin[idx];
                gbptr[idx].im =
                    fptr_z[idx].im * fftw_phase_cos[idx] -
                    fptr_z[idx].re * fftw_phase_sin[idx];
            }


            /*Do the backwards transform */
            fftwnd_one (p2, gbptr, beptr);
            /*This takes and stores the part of beta that is useful for this PE */
            assign_weight (sp, beptr, &beta_z[prjcount]);



            fptr += coarse_size;
            fptr_x += coarse_size;
            fptr_y += coarse_size;
            fptr_z += coarse_size;
            prjcount += ct.max_nlpoints;




        }                       /* end for ip */

        fftwnd_destroy_plan(p2);

    }                           /* end for ion */

    for (isp = 0; isp < ct.num_species; isp++)
    {
        fftw_free(ct.sp[isp].backward_wisdom);
    }

    my_free(beptr);
    my_free(fftw_phase_sin);

#if	DEBUG
    printf("PE: %d leave  get_nlop ...\n", pct.gridpe);
    fflush(NULL);
#endif

    if (pct.gridpe == 0)
    {

        printf(" get_nlop.c  done\n");

    }                           /* end if */
    /*    my_barrier(); */
    fflush(NULL);

}                               /* end get_nlop */
