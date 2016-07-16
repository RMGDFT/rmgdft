/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/



#include "portability.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include "Atomic.h"


void InitQfunct (std::unordered_map<std::string, InputKey *>& ControlMap)
{
    int idx, i, j, k, num, il, jl, ll;
    double *work = new double[MAX_RGRID]();
    double *qnmlig_tpr, *qnm_tpr;
    SPECIES *sp;
    char newname1[MAX_PATH];
    FILE *fqq = NULL;
    
    if(ct.norm_conserving_pp) return;

    Atomic *A = new Atomic();
    double *rgrid = A->GetRgrid();


    double *workr = new double[MAX_LOGGRID];

    for (int isp = 0; isp < ct.num_species; isp++)
    {

        sp = &ct.sp[isp];
        if (Verify ("write_pseudopotential_plots", true, ControlMap))
        {
            snprintf (newname1, MAX_PATH, "q_%s.xmgr", sp->atomic_symbol);
            if (pct.gridpe == 0)
            {
                if(NULL == (fqq = fopen (newname1, "w+")))
                    throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

            }
        }


        // Get range of q-functions. Assumes that all of them have roughly the same range
        qnm_tpr = sp->qnm;
        for (k = 0; k < sp->rg_points; k++)
        {
            if (sp->r[k] >= sp->rinner[0])
                work[k] = qnm_tpr[k];
            else
                work[k] = get_QnmL (0, 0, sp->r[k], sp);
        }
        sp->qradius = 2.5 * A->GetRange(work, sp->r, sp->rab, sp->rg_points);

        // Make adjustments so radii terminates on a grid point
        sp->qdim = Radius2grid (sp->qradius, ct.hmingrid/(double)Rmg_G->default_FG_RATIO);
        sp->qdim = sp->qdim/2*2 + 1;

        sp->qradius = 0.5 * ct.hmingrid * (double)(sp->qdim-1) / (double)Rmg_G->default_FG_RATIO;
        sp->qcut = 0.66 * sp->qradius;

        if (sp->qdim >= get_FNX_GRID()) sp->qdim = get_FNX_GRID();
        if (sp->qdim >= get_FNY_GRID()) sp->qdim = get_FNY_GRID();
        if (sp->qdim >= get_FNZ_GRID()) sp->qdim = get_FNZ_GRID();

        if (ct.max_Qpoints < (sp->qdim * sp->qdim * sp->qdim))
            ct.max_Qpoints = sp->qdim * sp->qdim * sp->qdim;

        num = (sp->nbeta * (sp->nbeta + 1)) * sp->nlc / 2;
        sp->qnmlig = new double[num * MAX_LOGGRID];
        idx = 0;

        for (i = 0; i < sp->nbeta; i++)
        {
            il = sp->llbeta[i];
            for (j = i; j < sp->nbeta; j++)
            {
                jl = sp->llbeta[j];
                idx = j * (j + 1) / 2 + i;
                qnm_tpr = sp->qnm + idx * MAX_RGRID;
                for (ll = abs (il - jl); ll <= il + jl; ll = ll + 2)
                {
                    for (k = 0; k < MAX_RGRID; k++)
                        work[k] = 0.0;

                    for (k = 0; k < sp->rg_points; k++)
                    {
                        if (sp->r[k] >= sp->rinner[ll])
                            work[k] = qnm_tpr[k]/sp->r[k]/sp->r[k];
                        else
                            work[k] = get_QnmL (idx, ll, sp->r[k], sp);
                    }
                    qnmlig_tpr = sp->qnmlig + (idx * sp->nlc + ll) * MAX_LOGGRID;

                    if (pct.gridpe == 0 && Verify ("write_pseudopotential_plots", true, ControlMap))
                    {
                        for (k = 0; k < sp->kkbeta; k++)
                            fprintf (fqq, "%e  %e\n", sp->r[k], work[k]);
                        fprintf (fqq, "&\n");

                    }

                    A->FilterPotential(work, sp->r, sp->rg_points, sp->qradius, ct.rhocparm, qnmlig_tpr,
                                        sp->rab, ll, sp->gwidth, sp->qcut, 40.0, ct.hmingrid/(double)Rmg_G->default_FG_RATIO);


                    /*Write final filtered Q function if requested*/
                    if (pct.gridpe == 0 && Verify ("write_pseudopotential_plots", true, ControlMap))
                    {

                        for (k = 0; k < MAX_LOGGRID; k++)
                        {
                            fprintf (fqq, "%e  %e\n", rgrid[k], qnmlig_tpr[k]);
                        }

                        fprintf (fqq, "&\n");
                    }
                }               /*end for ll */
            }                   /*end for j */
        }                       /*end for i */

        if (pct.gridpe == 0 && Verify ("write_pseudopotential_plots", true, ControlMap))
        {
            fclose (fqq);
        }

    }                           /*end for isp */ 

    delete A;
    delete [] workr;
    delete [] work;

}
