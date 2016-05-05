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


void InitQfunct (std::unordered_map<std::string, InputKey *>& ControlMap)
{
    int idx, i, j, k, num, il, jl, ll, it1, ibrav;
    double t1, t2;
    double work[MAX_RGRID];
    double *qnmlig_tpr, *drqnmlig_tpr, *qnm_tpr;
    SPECIES *sp;
    char newname1[MAX_PATH], newname2[MAX_PATH];
    FILE *fqq = NULL;
    FILE *fdq = NULL;
    
    if(ct.norm_conserving_pp) return;
    ibrav = Rmg_L.get_ibrav_type();

    Atomic *A = new Atomic();
    double *rgrid = A->GetRgrid();


    double *workr = new double[MAX_LOGGRID];

    double scale = 1.0;
    if (ibrav == CUBIC_BC)
        scale = 1.1;
    if (ibrav == CUBIC_FC)
        scale = 1.3;
    for (int isp = 0; isp < ct.num_species; isp++)
    {

        sp = &ct.sp[isp];
        if (Verify ("write_pseudopotential_plots", true, ControlMap))
        {
            snprintf (newname1, MAX_PATH, "q_%s.xmgr", sp->atomic_symbol);
            snprintf (newname2, MAX_PATH, "drq_%s.xmgr", sp->atomic_symbol);
            if (pct.gridpe == 0)
            {
                if(NULL == (fqq = fopen (newname1, "w+")))
                    throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

                if(NULL == (fdq = fopen (newname2, "w+")))
                    throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
            }
        }


        t1 = 2.0 * scale * (double) get_FG_RATIO() *sp->qradius / ct.hmingrid;
        /*		t1=2.0 *  scale * sp->qradius / ct.hmingrid;*/
        t1 = modf (t1, &t2);
        it1 = (int) t2;
        if (t1 > 0.5)
            it1++;
        if (!(it1 % 2))
            it1++;
        sp->qdim = it1;
        /*		sp->qdim = 2 * FG_NX * (it1 / 2) + 1;*/

        //        if ((sp->qdim >= get_FNX_GRID()) || (sp->qdim >= get_FNY_GRID())
        //            || (sp->qdim >= get_FNZ_GRID()))
        //            error_handler ("nlocal potential radius exceeds global grid size");
        if (sp->qdim >= get_FNX_GRID()) sp->qdim = get_FNX_GRID();
        if (sp->qdim >= get_FNY_GRID()) sp->qdim = get_FNY_GRID();
        if (sp->qdim >= get_FNZ_GRID()) sp->qdim = get_FNZ_GRID();

        if (ct.max_Qpoints < (sp->qdim * sp->qdim * sp->qdim))
            ct.max_Qpoints = sp->qdim * sp->qdim * sp->qdim;

        num = (sp->nbeta * (sp->nbeta + 1)) * sp->nlc / 2;
        sp->qnmlig = new double[num * MAX_LOGGRID];
        sp->drqnmlig = new double[num * MAX_LOGGRID];
        idx = 0;

        for (i = 0; i < sp->nbeta; i++)
        {
            il = sp->llbeta[i];
            for (j = i; j < sp->nbeta; j++)
            {
                jl = sp->llbeta[j];
                idx = j * (j + 1) / 2 + i;
                qnm_tpr = sp->qnm + idx * MAX_RGRID;
                for (ll = fabs (il - jl); ll <= il + jl; ll = ll + 2)
                {
                    for (k = 0; k < MAX_RGRID; k++)
                        work[k] = 0.0;

                    for (k = 0; k < sp->rg_points; k++)
                    {
                        if (sp->r[k] >= sp->rinner[ll])
                            work[k] = qnm_tpr[k];
                        else
                            work[k] = get_QnmL (idx, ll, sp->r[k], sp);
                    }
                    qnmlig_tpr = sp->qnmlig + (idx * sp->nlc + ll) * MAX_LOGGRID;
                    drqnmlig_tpr = sp->drqnmlig + (idx * sp->nlc + ll) * MAX_LOGGRID;

                    if (pct.gridpe == 0 && Verify ("write_pseudopotential_plots", true, ControlMap))
                    {
                        for (k = 0; k < sp->kkbeta; k++)
                            fprintf (fqq, "%e  %e\n", sp->r[k], work[k]);
                        fprintf (fqq, "&\n");

                    }

                    A->FilterPotential(work, sp->r, sp->rg_points, sp->qradius, 0.25, ct.qcparm, qnmlig_tpr,
                                        sp->rab, ll, sp->gwidth, sp->nlrcut[sp->llbeta[i]], sp->rwidth, drqnmlig_tpr, sp->qdim);

                    /*Is this necessary ???*/
                    if (ll)
                        qnmlig_tpr[0] = 0.0;


                    /*Write final filtered Q function if requested*/
                    if (pct.gridpe == 0 && Verify ("write_pseudopotential_plots", true, ControlMap))
                    {

                        for (k = 0; k < MAX_LOGGRID; k++)
                        {
                            fprintf (fqq, "%e  %e\n", rgrid[k], qnmlig_tpr[k]);
                            fprintf (fdq, "%e  %e\n", rgrid[k], drqnmlig_tpr[k]);
                        }

                        fprintf (fqq, "&\n");
                        fprintf (fdq, "&\n");
                    }
                }               /*end for ll */
            }                   /*end for j */
        }                       /*end for i */

        if (pct.gridpe == 0 && Verify ("write_pseudopotential_plots", true, ControlMap))
        {
            fclose (fqq);
            fclose (fdq);
        }

    }                           /*end for isp */ 

    delete A;
    delete [] workr;

}
