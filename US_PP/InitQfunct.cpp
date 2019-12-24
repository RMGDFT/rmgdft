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


void InitQfunct ()
{
    if(ct.norm_conserving_pp) return;
    int idx, i, j, k, num, il, jl, ll;
    double *qnmlig_tpr, *qnm_tpr;
    SPECIES *sp;
    char newname1[MAX_PATH];
    FILE *fqq = NULL;
    


    Atomic *A = new Atomic();
    double *rgrid = A->GetRgrid();



    for (int isp = 0; isp < ct.num_species; isp++)
    {

        sp = &Species[isp];
        if(!std::strcmp(sp->atomic_symbol, "DLO")) continue;
        if (ct.write_pp_flag)
        {
            snprintf (newname1, MAX_PATH, "q_%s.xmgr", sp->atomic_symbol);
            if (pct.gridpe == 0)
            {
                if(NULL == (fqq = fopen (newname1, "w+")))
                    throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

            }
        }

        //  change from qnm to qnm_l
        if(!sp->q_with_l)
        {
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

                        qnmlig_tpr = sp->qnmlig + (idx * sp->nlc + ll) * MAX_LOGGRID;

                        for (k = 0; k < MAX_RGRID; k++)
                            sp->qnm_l[k + (idx * sp->nlc + ll) * MAX_RGRID ] = 0.0;

                        for (k = 0; k < sp->rg_points; k++)
                        {
                            sp->qnm_l[k + (idx * sp->nlc + ll) * MAX_RGRID ] = qnm_tpr[k]/sp->r[k]/sp->r[k];
                            if(sp->nqf > 0 && sp->r[k] < sp->rinner[ll])
                                sp->qnm_l[k + (idx * sp->nlc + ll) * MAX_RGRID ] = get_QnmL (idx, ll, sp->r[k], sp);
                        }
                    }
                }
            }
        }
        else
        {
            for(i = 0; i < (sp->nlc * sp->nbeta *(sp->nbeta +1)) /2; i++)
            for (k = 0; k < sp->rg_points; k++)
            {
                sp->qnm_l[k + i * MAX_RGRID ] /= sp->r[k]*sp->r[k];
            }
        }

        double qradius = ct.min_qradius;
        for(i = 0; i < (sp->nlc * sp->nbeta *(sp->nbeta +1)) /2; i++)
        {
            qnm_tpr = &sp->qnm_l[i * MAX_RGRID];
            sp->qradius = 2.5 * A->GetRange(qnm_tpr, sp->r, sp->rab, sp->rg_points, 0.999999999);
            sp->qradius = std::max(sp->qradius, qradius);
            qradius = sp->qradius;

        }

        sp->qradius = std::min(sp->qradius, ct.max_qradius);

        // Make adjustments so radii terminates on a grid point
        sp->qdim = Radius2grid (sp->qradius, ct.hmingrid/(double)Rmg_G->default_FG_RATIO, Rmg_L.get_ibrav_type(), false);
        sp->qdim = sp->qdim/2*2 + 1;

        sp->qcut = 0.5 * sp->qradius;

        if (sp->qdim >= get_FNX_GRID()) sp->qdim = get_FNX_GRID();
        if (sp->qdim >= get_FNY_GRID()) sp->qdim = get_FNY_GRID();
        if (sp->qdim >= get_FNZ_GRID()) sp->qdim = get_FNZ_GRID();

        if (ct.max_Qpoints < (sp->qdim * sp->qdim * sp->qdim))
            ct.max_Qpoints = sp->qdim * sp->qdim * sp->qdim;

        num = (sp->nbeta * (sp->nbeta + 1)) * sp->nlc / 2;
        sp->qnmlig = new double[num * MAX_LOGGRID]();
        idx = 0;

        for (i = 0; i < sp->nbeta; i++)
        {
            il = sp->llbeta[i];
            for (j = i; j < sp->nbeta; j++)
            {
                jl = sp->llbeta[j];
                idx = j * (j + 1) / 2 + i;
                for (ll = abs (il - jl); ll <= il + jl; ll = ll + 2)
                {

                    qnmlig_tpr = sp->qnmlig + (idx * sp->nlc + ll) * MAX_LOGGRID;
                    qnm_tpr = sp->qnm_l + (idx * sp->nlc + ll) *  MAX_RGRID;

                    if (pct.gridpe == 0 && ct.write_pp_flag)
                    {
                        for (k = 0; k < sp->kkbeta; k++)
                            fprintf (fqq, "%e  %e\n", sp->r[k], qnm_tpr[k]);
                        fprintf (fqq, "&\n");
                    }

                    A->FilterPotential(qnm_tpr, sp->r, sp->rg_points, sp->qradius, ct.rhocparm, qnmlig_tpr,
                            sp->rab, ll, sp->gwidth, sp->qcut, 1.0, ct.hmingrid/(double)Rmg_G->default_FG_RATIO);

                    /*Write final filtered Q function if requested*/
                    if (pct.gridpe == 0 && ct.write_pp_flag)
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

        // Raw q function from pp is no longer needed so free it's memory
        delete []  sp->qnm;

        if (pct.gridpe == 0 && ct.write_pp_flag)
        {
            fclose (fqq);
        }

    }                           /*end for isp */ 

    delete A;

}
