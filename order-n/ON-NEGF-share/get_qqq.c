/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "md.h"

void get_qqq()
{
    int idx, i, j, ion, ispe;
    int nh, ncount, icount;
    REAL *qnmI, *qqq, sum;
    ION *iptr;
    SPECIES *sp;
    FILE *ftpr;

    if (pct.gridpe == 0)
        ftpr = fopen("check_q.txt", "w+");


    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        nh = sp->num_projectors;
        ncount = pct.Qidxptrlen[ion];
        qnmI = pct.augfunc[ion];

        if (pct.qqq[ion] == NULL)
            my_calloc( pct.qqq[ion], nh * nh, REAL );
        qqq = pct.qqq[ion];

        if (pct.gridpe == 0)
            fprintf(ftpr, "%% for ion %d :\n", ion);

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
                sum = real_sum_all(sum, pct.grid_comm);
                sum = sum * ct.vel_f;
                if (fabs(sum) < 1.0e-8)
                    sum = 0.0;
                if (pct.gridpe == 0)
                    fprintf(ftpr, "i=%d j=%d q_cal=%15.8f q_rel=%15.8f\n", i,
                            j, sum, sp->qqq[i][j]);
                qqq[i * nh + j] = sum;
                if (i != j)
                    qqq[j * nh + i] = qqq[i * nh + j];

                idx++;
            }                   /*end for j */
        }                       /*end for i */
        if (pct.gridpe == 0)
            fprintf(ftpr, "\n");
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
    if (pct.gridpe == 0)
        fclose(ftpr);
}
