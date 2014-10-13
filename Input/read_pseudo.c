/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/read_pseudo.c *****
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
 *   void read_pseudo(void)
 *   Reads in species and pseudopotential information.
 * INPUTS
 *   no explicit input
 * OUTPUT
 *   pseudopotential informations are stored in ct.sp...
 * PARENTS
 *   main.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "main.h"


void read_pseudo (void)
{

    int i, j, k, idx, indx, idx1, idx2, idx3, idx4, nmb;
    int l, ih, check;
    rmg_double_t  ddd0[6][6], qqq[6][6];
    /*rmg_double_t ddd[6][6]; */
    int max_nlprojectors = 0, nlc;
    SPECIES *sp;
    char tbuf[MAX_CHAR], *tptr, *endptr;

    ct.max_l = 0;


    for (i = 0; i < ct.num_species; i++)
    {

        sp = &ct.sp[i];

        // Allocate arrays
        my_malloc (sp->r, MAX_RGRID, double);
        my_malloc (sp->rab, MAX_RGRID, double);
        my_malloc (sp->vloc0, MAX_RGRID, double);
        my_malloc (sp->cr, MAX_RGRID, double);
        my_malloc (sp->rspsco, MAX_RGRID, double);

        get_data (sp->pseudo_filename, NULL, LINES|INIT, NULL);

//        my_malloc(sp->description, MAX_CHAR, char);

        /* Description */
        get_data (sp->pseudo_filename, &sp->description, ITEM | STR, "ultrasoft pseudo-potential"); 
        //dprintf( "desc=:%s:\n", sp->description );


        /* Read in the nlcc flag */
        get_data (sp->pseudo_filename, &sp->nlccflag, ITEM | INT, NULL); 

        /* Read in the atomic number */
        get_data (sp->pseudo_filename, &sp->atomic_number, ITEM | INT, NULL); 
        /*Determine atomic symbol */
        sp->atomic_symbol = get_symbol (sp->atomic_number);

        /* Read in the atomic mass */
        get_data (sp->pseudo_filename, &sp->atomic_mass, ITEM | DBL, NULL); 

        /* Read in the number of valence electrons */
        get_data (sp->pseudo_filename, &sp->zvalence, ITEM | DBL, NULL); 

        /* Gaussian charge parameter */
        get_data (sp->pseudo_filename, &sp->rc, ITEM | DBL, NULL); 

        /* Number of potentials */
        get_data (sp->pseudo_filename, &sp->num_potentials, ITEM | INT, NULL); 

        if (sp->num_potentials > (MAX_L + 1))
            error_handler ("%d is more than (MAX_L + 1)=%d potentials", sp->num_potentials, MAX_L + 1);


        /* l-values for the potentials */
        for (j = 0; j < sp->num_potentials; j++)
        {

            get_data (sp->pseudo_filename, &sp->lval[j], ITEM | INT, NULL); 
            if (sp->lval[j] > MAX_L)
                error_handler ("Maximum l-value too large");
        }

        /* L-value for the local potential */
        get_data (sp->pseudo_filename, &sp->local, ITEM | INT, NULL); 
        if (sp->local > MAX_L)
            error_handler ("Maximum l-value too large");

        /* Local potential radius */
        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL);
        idx = sscanf ( tbuf, " %lf %lf %lf %lf %lf ", &sp->lradius, &sp->aradius, &sp->acut, &sp->agwidth, &sp->arwidth ); 
        if (idx == 1)
        {
            sp->aradius = 9.0;
            sp->acut = 7.0;
            sp->agwidth = 10.0;
            sp->arwidth = 25.0;
        }

        if ((idx != 1) && (idx != 5))
            error_handler("Unexpected number of parameters (%d), only 1 or 5 are valid. The string was %s", idx, tbuf);

        /* Non-Local potential radius */
        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
        idx = sscanf ( tbuf, " %lf %lf ", &sp->nlradius, &sp->qradius ); 
        if ( 2 != idx )
            error_handler ( "Should read in 2 doubles, got %d here, nlradius and qradius", idx );

        /* Number of points in the radial grid */
        get_data (sp->pseudo_filename, &sp->rg_points, ITEM | INT, NULL); 

        /*Number of beta function and the number of point in radial grid */
        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
        if ( 2 != sscanf ( tbuf, " %d %d ", &sp->nbeta, &sp->kkbeta ) )
            error_handler ( "Should read in 2 ints here, nbeta and kkbeta");
        if (sp->nbeta > MAX_NB)
            error_handler ("too may beta functions");

        /* Filtering information */
        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
        sp->lrcut = strtod( tbuf, &tptr);
        for (j = 0; j < sp->num_potentials; j++)
        {
            if (sp->lval[j] != sp->local)
            {
                sp->nlrcut[sp->lval[j]] = strtod( tptr, &endptr);
                if ( endptr == tptr )
                    error_handler ( "Missing nonlocal cutoff radius information");
                else 
                    tptr = endptr;
            }
        }
        
        if ( 2 != sscanf ( tptr, " %lf %lf ", &sp->rwidth, &sp->gwidth ) )
                error_handler ( "Should read in 2 floats here, sp->rwidth, sp->gwidth ");
        Dprintf( "sp->rwidth=%lf, sp->gwidth=%lf", sp->rwidth, sp->gwidth);

        /*read in the nlc and nqf */
        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
        if ( 2 != sscanf ( tbuf, " %d %d ", &sp->nlc, &sp->nqf  ) )
            error_handler ( "Should read in 2 ints here, nlc and nqf");

        /*read in the rinner for Q_L(r) function */
        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
        tptr = tbuf;
         my_malloc(sp->rinner, sp->nlc, int);
        for (j = 0; j < sp->nlc; j++)
        {
            sp->rinner[j] = strtod( tptr, &endptr);
            if ( endptr == tptr )
                error_handler ( "Missing inner cutoff radius information");
            else 
                tptr = endptr;
        }
        /*read in the Matrix ddd0(nbeta,nbeta) */
        for (j = 0; j < sp->nbeta; j++)
        {
            get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
            tptr = tbuf;
            for (k = 0; k < sp->nbeta; k++)
            {
                ddd0[j][k] = strtod( tptr, &endptr);
                if ( endptr == tptr )
                    error_handler ( "Missing ddd matrix information");
                else 
                    tptr = endptr;
            }
        }

        /*read in the matrix qqq(nbeta,nbeta) */
        for (j = 0; j < sp->nbeta; j++)
        {
            get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
            tptr = tbuf;
            for (k = 0; k < sp->nbeta; k++)
            {
                qqq[j][k] = strtod( tptr, &endptr);
                if ( endptr == tptr )
                    error_handler ( "Missing qqq matrix information");
                else 
                    tptr = endptr;
            }
        }

        Dprintf( "Coefficients of the pseudoized Q_L(r)  Start\n", tbuf);
        /*read in the coefficient of the pseudoized Q_L(r) function */
        nlc = (sp->nbeta * (sp->nbeta + 1)) / 2;
        my_malloc (sp->qfcoef, nlc * sp->nlc * sp->nqf, rmg_double_t);
        idx = 0;
        for (idx1 = 0; idx1 < sp->nbeta; idx1++)
        {
            for (idx2 = idx1; idx2 < sp->nbeta; idx2++)
            {
                nmb = idx2 * (idx2 + 1) / 2 + idx1;
                for (idx3 = 0; idx3 < sp->nlc; idx3++)
                {
                    for (idx4 = 0; idx4 < sp->nqf; idx4++)
                    {
                        indx = nmb * sp->nlc * sp->nqf + idx3 * sp->nqf + idx4;
                        if ( idx%4 == 0 )
                        {
                            get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
                            tptr = tbuf;
                        }
                        Dprintf( "casting :%s: into pseudized Q_L(r) at index %d/%d\n", tbuf, indx, nlc * sp->nlc * sp->nqf);
                        sp->qfcoef[indx] = strtod ( tptr, &endptr );
                        idx++;
                        tptr = endptr;
                    }
                }
            }
        }


        /* Next read in the radial grid,rab(r), Vloc(r), and Vloc_ion(r) for each species */
        Dprintf( "Radial Grid Start\n", tbuf);
        for (j = 0; j < sp->rg_points; j++)
        {

            get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
            Dprintf( "Radial Grid index is %d/%d\n", j, sp->rg_points);
            if ( 3 != sscanf ( tbuf, " %lf %lf %lf ", &sp->r[j], &sp->rab[j], &sp->vloc0[j]  ) )
                error_handler ( "Should read in 3 floats here, r, rab, vloc0");
        }

        sp->aa = (sp->r[0] * sp->r[0]) / (sp->r[1] - 2 * sp->r[0]);
        sp->bb = log (sp->r[1] / sp->r[0] - 1);

        for (j = 0; j < sp->nbeta; j++)
        {
            get_data (sp->pseudo_filename, &sp->llbeta[j], ITEM | INT, NULL); 
            
            /*Make sure that beta function has nlrcut defined*/ 
            check = 0;
            for (k = 0; k < sp->num_potentials; k++)
            {
                if (sp->lval[k] != sp->local)
                {
                    if (sp->lval[k] == sp->llbeta[j])
                        check = 1;
                }
            }

            if (!check)
                error_handler("Beta function with l=%d has undefined nlrcut", sp->llbeta[j]);
            if (sp->llbeta[j] > ct.max_l)
                ct.max_l = sp->llbeta[j];
            // Allocate memory for the beta
            my_malloc(sp->beta[j], MAX_RGRID, double);


            for (k = 0; k < MAX_RGRID; k++)
                sp->beta[j][k] = 0.0;
            for ( k = 0; k < sp->kkbeta; k+=4 )
            {
                get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
                idx = sscanf ( tbuf, " %lf %lf %lf %lf ",
                        &sp->beta[j][k + 0],
                        &sp->beta[j][k + 1],
                        &sp->beta[j][k + 2],
                        &sp->beta[j][k + 3]); 
#if DEBUG
                if ( idx != 4 )
                { /* Should normally read in 4 doubles per line */
                    if ( sp->kkbeta - k < 4 )
                    { /* except, perhaps the last line */
                        if ( idx != sp->kkbeta%4 )
                        { /* Should only read in the remainder of kkbeta%4 */
                            error_handler ( "Should read in %d beta floats here", sp->kkbeta%4 );
                        }
                    } else { /* Not the last line and did not get 4 values */
                        error_handler ( "Should have read in 4 beta floats here" );
                    }
                }
#endif
            }

            /* Beta function is defined on only a subset of sp->rg_points, but filtering uses all sp->rg_points, so 
             * we initialize the remaining data to zero */
            for ( k = sp->kkbeta; k < sp->rg_points; k++)
                sp->beta[j][k] = 0.0;
        }
        my_calloc (sp->qnm, nlc * MAX_RGRID, rmg_double_t);
        for (idx = 0; idx < nlc; idx++)
        {
            for (idx1 = 0; idx1 < MAX_RGRID; idx1++)
                sp->qnm[idx * MAX_RGRID + idx1] = 0.0;
        }
        if(sp->nqf != 0) {
            for (idx1 = 0; idx1 < sp->nbeta; idx1++)
            {
                for (idx2 = idx1; idx2 < sp->nbeta; idx2++)
                {
                    nmb = idx2 * (idx2 + 1) / 2 + idx1;
                    for ( idx3 = 0; idx3 < sp->kkbeta; idx3+=4 )
                    {
                        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
                        idx = sscanf ( tbuf, " %lf %lf %lf %lf ",
                                &sp->qnm[nmb * MAX_RGRID + idx3 + 0],
                                &sp->qnm[nmb * MAX_RGRID + idx3 + 1],
                                &sp->qnm[nmb * MAX_RGRID + idx3 + 2],
                                &sp->qnm[nmb * MAX_RGRID + idx3 + 3]);
                    }
                }
            }
        }

        if (sp->nlccflag)
        {

            /* Read in the pseudo valence and core charge densities */
            for (k = 0; k < sp->rg_points; k+=4)
            {
                get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
                idx = sscanf ( tbuf, " %lf %lf %lf %lf ",
                        &sp->rspsco[k + 0], 
                        &sp->rspsco[k + 1], 
                        &sp->rspsco[k + 2], 
                        &sp->rspsco[k + 3]); 
            }
        }

        for (j = 0; j < 18; j++)
        {
            sp->nhtol[j] = 0;
            sp->nhtom[j] = 0;
            sp->indv[j] = 0;
        }

        ih = 0;
        for (j = 0; j < sp->nbeta; j++)
        {
            l = sp->llbeta[j];
            for (k = 0; k < 2 * l + 1; k++)
            {
                sp->nhtol[ih] = l;
                sp->nhtom[ih] = k;
                sp->indv[ih] = j;
                ++ih;
            }
        }
        sp->nh = ih;
        if (ih > max_nlprojectors)
            max_nlprojectors = ih;
        if (max_nlprojectors > MAX_NL)
            error_handler ("too many nonlocal projectors");


        for (j = 0; j < ih; j++)
        {
            for (k = 0; k < ih; k++)
            {
                if ((sp->nhtol[j] == sp->nhtol[k]) && (sp->nhtom[j] == sp->nhtom[k]))
                {
                    sp->ddd0[j][k] = ddd0[sp->indv[j]][sp->indv[k]];
                    sp->qqq[j][k] = qqq[sp->indv[j]][sp->indv[k]];
                    /*                          sp->ddd[j][k]=ddd[indv[j]][indv[k]]; */
                }
                else
                {
                    sp->ddd0[j][k] = 0.0;
                    sp->qqq[j][k] = 0.0;
                    /*                          sp->ddd[j][k]=0.0; */
                }
            }
        }


        /*Read atomic wavefunctions */
        if (verify ("start_mode","LCAO Start"))
        {


            /*Milliken radius */
            sp->mill_radius = 9.0;

            sp->num_atomic_waves = 0;

            /* Number of wavefunctions (there should be one for s, p, d etc. states) 
             * We do conditional reading so that ON LCAO start, which does not use wavefunction from PP file
             * can go on */
            if (get_data (sp->pseudo_filename, &sp->num_atomic_waves, ITEM | INT, NULL)) 
            {
                
                my_malloc(sp->atomic_wave, sp->num_atomic_waves, rmg_double_t *);
                my_malloc(sp->awave_lig, sp->num_atomic_waves, rmg_double_t *);

                for (j = 0; j < sp->num_atomic_waves; j++ )
                {
                    /*Allocate and zero memory for atomic wave functions */
                    my_malloc(sp->atomic_wave[j], sp->rg_points, rmg_double_t);
                    my_malloc(sp->awave_lig[j], MAX_LOCAL_LIG, rmg_double_t);

                    get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 

                    idx = sscanf (tbuf, " %s %d %lf",&sp->atomic_wave_label[j][0], &sp->atomic_wave_l[j], &sp->atomic_wave_oc[j]);

                    /* IS THIS NECESSARY ?*/
                    /*for (k = 0; k < sp->rg_points; k++)
                      sp->atomic_wave[j][k] = 0.0;*/

                    for (k = 0; k < sp->rg_points; k+=4)
                    {
                        get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
                        idx = sscanf ( tbuf, " %lf %lf %lf %lf ",
                                &sp->atomic_wave[j][k + 0], 
                                &sp->atomic_wave[j][k + 1], 
                                &sp->atomic_wave[j][k + 2], 
                                &sp->atomic_wave[j][k + 3]); 
                    }
                }


                /*Read atomic charge density*/
                my_malloc(sp->atomic_rho, sp->rg_points, rmg_double_t);
                for (k = 0; k < sp->rg_points; k+=4)
                {
                    get_data (sp->pseudo_filename, tbuf, ITEM | STR, NULL); 
                    idx = sscanf ( tbuf, " %lf %lf %lf %lf ",
                            &sp->atomic_rho[k + 0], 
                            &sp->atomic_rho[k + 1], 
                            &sp->atomic_rho[k + 2], 
                            &sp->atomic_rho[k + 3]); 
                }

            }

        }                       /*end if (ct.domilliken) */


        /* Remove data node, it is not needed anymore. */
        get_data (sp->pseudo_filename, NULL, END|INIT, NULL);
    }                           /* end for(i = 0;i < ct.num_species;i++) */


    /* Set the maximum number of non-local projecters needed */
    ct.max_nl = max_nlprojectors;

}                               /* end read_pseudo */



/******/
