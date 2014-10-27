/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"

void get_qqq ()
{
    int idx, i, j, ion, ispe;
    int nh, ncount, icount;
    rmg_double_t *qnmI, *qqq, sum;
    ION *iptr;
    SPECIES *sp;
    FILE *ftpr;
    char filename[MAX_PATH];
    const bool SET=true;

    if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
    {
	snprintf (filename, MAX_PATH, "%s.q.txt", ct.basename);
        my_fopen (ftpr, filename, "w+");
    }


    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        nh = sp->nh;
        ncount = pct.Qidxptrlen[ion];
        qnmI = pct.augfunc[ion];

        if (pct.qqq[ion] == NULL)
            my_malloc (pct.qqq[ion], nh * nh, rmg_double_t);
        qqq = pct.qqq[ion];

        if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
            fprintf (ftpr, "%% for ion %d :\n", ion);

        idx = 0;
        for (i = 0; i < nh; i++)
        {
            for (j = i; j < nh; j++)
            {
                sum = 0.0;
                if (ncount)
                {
                    qnmI = pct.augfunc[ion] + idx * ncount;
                    for (icount = 0; icount < ncount; icount++)
                    {
                        sum += qnmI[icount];
                    }
                }
                sum = real_sum_all (sum, pct.grid_comm);
                sum = sum * get_vel_f();
                if (fabs (sum) < 1.0e-8)
                    sum = 0.0;
                if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
                    fprintf (ftpr, "i=%d j=%d q_cal=%15.8f q_rel=%15.8f\n",
                             i, j, sum, sp->qqq[i][j]);

                qqq[i * nh + j] = sum;
                if (i != j)
                    qqq[j * nh + i] = qqq[i * nh + j];

                idx++;
            }                   /*end for j */
        }                       /*end for i */
        if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
            fprintf (ftpr, "\n");
    }                           /*end for ion */
/*
	for(i=0;i<ct.num_species;i++) {
		if(pct.gridpe==0) {
			fprintf(ftpr,"%% for Specie %d :\",i);
			fprintf(ftpr,"N    M    ");
			for(j=0;j<ct.num_ions,j++) 
				iptr=&ct.ions[ion];
				if(i==iptr->species) fprintf(ftpr,"ION%d   ",j);
			}
	*/
    if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
        fclose (ftpr);
}
