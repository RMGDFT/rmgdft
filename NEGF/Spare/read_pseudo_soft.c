/************************** SVN Revision Information **************************
 **    $Id: read_pseudo_soft.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/read_softpseudo_soft.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 1995  Shuchun Wang,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void read_softpseudo(void)
 *   Reads in species and pseudopotential information.
 * INPUTS
 *   no explicit input
 * OUTPUT
 *   pseudopotential informations are stored in ct.sp...
 * PARENTS
 *   md.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "md.h"





char *get_num (char *str);
char *get_line (char *buf, FILE * fh);


void read_pseudo_soft (void)
{

    int i, j, k, idx, indx, idx1, idx2, idx3, idx4, nmb;
    int l, ih, indv[18], nhtol[18], nhtom[18];
    REAL temp, ddd0[6][6], qqq[6][6], ddd[6][6];
    int max_num_projectors = 0, nlc;
    FILE *fhand;
    SPECIES *sp;
    char tbuf[200], *tptr;

    if (NULL == (fhand = fopen (ct.pspfile, "rw")))
        error_handler ("Unable to open pseudopotential file");


    for (i = 0; i < ct.num_species; i++)
    {


        sp = ct.sp[i];

        /* Description */
        strcpy (sp->description, get_line (tbuf, fhand));


        /* Read in the nlcc flag */
        sp->nlccflag = atoi (get_line (tbuf, fhand));

        /* Read in the atomic number */
        sp->atomic_number = atoi (get_line (tbuf, fhand));

        /*Determine atomic symbol */
        sp->atomic_symbol = get_symbol (sp->atomic_number);

        /* Read in the atomic mass */
        sp->atomic_mass = atof (get_line (tbuf, fhand));

        /* Read in the number of valence electrons */
        sp->zvalence = atof (get_line (tbuf, fhand));

        /* Gaussian charge parameter */
        sp->rc = atof (get_line (tbuf, fhand));

        /* Number of potentials */
        sp->num_potentials = atoi (get_line (tbuf, fhand));

        if (sp->num_potentials > (MAX_L + 1))
            error_handler ("Too many potentials");


        /* l-values for the potentials */
        for (j = 0; j < sp->num_potentials; j++)
        {

            sp->lval[j] = atoi (get_line (tbuf, fhand));
            if (sp->lval[j] > MAX_L)
                error_handler ("Maximum l-value too large");


        }                       /* end for */

        /* L-value for the local potential */
        sp->local = atoi (get_line (tbuf, fhand));
        if (sp->local > MAX_L)
            error_handler ("Maximum l-value too large");

        /* Local potential radius */
        sp->lradius = atof (get_line (tbuf, fhand));

        /* Nonlocal potential radius */
        sp->nlradius = atof (get_line (tbuf, fhand));
        sp->qradius = atof (get_num (tbuf));

        /* Number of points in the radial grid */
        sp->rg_points = atoi (get_line (tbuf, fhand));

        /*Number of beta function and the number of point in radial grid */
        sp->nbeta = atoi (get_line (tbuf, fhand));
        sp->kkbeta = atoi (get_num (tbuf));
        if (sp->nbeta > MAX_NB)
            error_handler ("too may beta function");

        /* Filtering information */
        tptr = get_line (tbuf, fhand);
        sp->lrcut = atof (tptr);
        for (j = 0; j < sp->num_potentials; j++)
        {
            if (sp->lval[j] != sp->local)
            {
                if (NULL == (tptr = get_num (tptr)))
                    error_handler ("Missing nonlocal cutoff radius information");
                sp->nlrcut[sp->lval[j]] = atof (tptr);
            }
        }
        if (NULL == (tptr = get_num (tptr)))
            error_handler ("Missing filtering information");
        sp->rwidth = atof (tptr);
        if (NULL == (tptr = get_num (tptr)))
            error_handler ("Missing filtering information");
        sp->gwidth = atof (tptr);

        /*read in the nlc and nqf */
        sp->nlc = atoi (get_line (tbuf, fhand));
        sp->nqf = atoi (get_num (tbuf));

        /*read in the rinner for Q_L(r) function */
        tptr = get_line (tbuf, fhand);
        for (j = 0; j < sp->nlc; j++)
        {
            sp->rinner[j] = atof (tptr);
            tptr = get_num (tptr);
        }

        /*read in the Matrix ddd0(nbeta,nbeta) */
        for (j = 0; j < sp->nbeta; j++)
        {
            tptr = get_line (tbuf, fhand);
            for (k = 0; k < sp->nbeta; k++)
            {
                ddd0[j][k] = atof (tptr);
                tptr = get_num (tptr);
            }
        }

        /*read in the Matrix ddd(nbeta,nbeta) */
/*      for(j=0;j<sp->nbeta; j++){
	      tptr=get_line(tbuf,fhand);
	      for(k=0;k<sp->nbeta;k++) {
		      ddd[j][k]=atof(tptr);
		      tptr=get_num(tptr);
	      }
      }
*/

        /*read in the matrix qqq(nbeta,nbeta) */
        for (j = 0; j < sp->nbeta; j++)
        {
            tptr = get_line (tbuf, fhand);
            for (k = 0; k < sp->nbeta; k++)
            {
                qqq[j][k] = atof (tptr);
                tptr = get_num (tptr);
            }
        }

        /*read in the coefficient of the pseudoized Q_L(r) function */
        nlc = (sp->nbeta * (sp->nbeta + 1)) / 2;
        my_malloc_init( sp->qfcoef, nlc * sp->nlc * sp->nqf, REAL );
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

                        if ((idx) % 4 == 0)
                            tptr = get_line (tbuf, fhand);
                        temp = atof (tptr);
                        sp->qfcoef[indx] = temp;
                        idx++;
                        tptr = get_num (tptr);
                    }
                }
            }
        }


        /* Next read in the radial grid,rab(r), Vloc(r), and Vloc_ion(r) for each specie */
        for (j = 0; j < sp->rg_points; j++)
        {

            sp->r[j] = atof (get_line (tbuf, fhand));

            if (NULL == (tptr = get_num (tbuf)))
                error_handler ("Bad input field in pseudopotential file");
            sp->rab[j] = atof (tptr);


            if (NULL == (tptr = get_num (tptr)))
                error_handler ("Bad input field in pseudopotential file");

            sp->vloc0[j] = atof (tptr);

        }                       /* end for */

        sp->aa = (sp->r[0] * sp->r[0]) / (sp->r[1] - 2 * sp->r[0]);
        sp->bb = log (sp->r[1] / sp->r[0] - 1);

        for (j = 0; j < sp->nbeta; j++)
        {
            sp->llbeta[j] = atoi (get_line (tbuf, fhand));
            for (k = 0; k < MAX_RGRID; k++)
                sp->beta[j][k] = 0.0;
            for (k = 0; k < (sp->kkbeta) / 4; k++)
            {
                idx = k * 4;
                sp->beta[j][idx + 0] = atof (get_line (tbuf, fhand));
                sp->beta[j][idx + 1] = atof (tptr = get_num (tbuf));
                sp->beta[j][idx + 2] = atof (tptr = get_num (tptr));
                sp->beta[j][idx + 3] = atof (get_num (tptr));

            }
            if (((sp->kkbeta) % 4) != 0)
            {
                idx = k * 4;
                tptr = get_line (tbuf, fhand);
                for (k = 0; k < (sp->kkbeta) % 4; k++)
                {
                    sp->beta[j][idx + k] = atof (tptr);
                    tptr = get_num (tptr);
                }
            }
        }

        my_malloc_init( sp->qnm, nlc * MAX_RGRID, REAL );
        for (idx = 0; idx < nlc; idx++)
        {
            for (idx1 = 0; idx1 < MAX_RGRID; idx1++)
                sp->qnm[idx * MAX_RGRID + idx1] = 0.0;
        }
        for (idx1 = 0; idx1 < sp->nbeta; idx1++)
        {
            for (idx2 = idx1; idx2 < sp->nbeta; idx2++)
            {
                nmb = idx2 * (idx2 + 1) / 2 + idx1;
                for (idx3 = 0; idx3 < sp->kkbeta; idx3++)
                {
                    if ((idx3 % 4) == 0)
                        tptr = get_line (tbuf, fhand);
                    sp->qnm[nmb * MAX_RGRID + idx3] = atof (tptr);
                    /*                    printf("qnm[%d][%d]=%e\n",nmb,idx3,sp->qnm[nmb][idx3]); */
                    tptr = get_num (tptr);
                }
            }
        }


        if (sp->nlccflag)
        {

            /* Read in the pseudo valence and core charge densities */
            for (k = 0; k < sp->rg_points; k++)
            {

                /* Radial grid for the densities */
                if ((k % 4) == 0)
                    tptr = get_line (tbuf, fhand);
                sp->rspsco[k] = atof (tptr);
                tptr = get_num (tptr);

            }                   /* end for */

        }                       /* end if */


        /* Radial grid for the densities */
/*      for(k = 0;k < sp->rg_points;k++) {

           if((k%4)==0) tptr=get_line(tbuf,fhand);
           sp->rsatom[k] = atof(tptr);
           tptr=get_num(tptr);

       } 
*/

        ih = 0;
        for (j = 0; j < sp->nbeta; j++)
        {
            l = sp->llbeta[j];
            for (k = 0; k < 2 * l + 1; k++)
            {
                nhtol[ih] = l;
                nhtom[ih] = k;
                indv[ih] = j;
                ++ih;
            }
        }

        sp->num_projectors = ih;
        if (ih > max_num_projectors)
            max_num_projectors = ih;
        if (max_num_projectors > MAX_NL)
            error_handler ("too many nonlocal projectors");


        for (j = 0; j < ih; j++)
        {
            for (k = 0; k < ih; k++)
            {
                if ((nhtol[j] == nhtol[k]) && (nhtom[j] == nhtom[k]))
                {
                    sp->ddd0[j][k] = ddd0[indv[j]][indv[k]];
                    sp->qqq[j][k] = qqq[indv[j]][indv[k]];
                    /*                    sp->ddd[j][k]=ddd[indv[j]][indv[k]]; */
                }
                else
                {
                    sp->ddd0[j][k] = 0.0;
                    sp->qqq[j][k] = 0.0;
                    /*     sp->ddd[j][k]=0.0; */
                }
            }
        }

    }                           /* end for */

    fclose (fhand);

    /* Set the maximum number of non-local projecters needed */
    ct.max_nl = max_num_projectors;
}                               /* end read_species */

/******/
