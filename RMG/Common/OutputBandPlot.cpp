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
#include <stdlib.h>
#include <stdio.h>
#include "portability.h"
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "RmgException.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "../Headers/prototypes.h"
#include "Plots.h"


template void OutputBandPlot(Kpoint<double> **);
template void OutputBandPlot(Kpoint<std::complex<double> > **);

template <typename KpointType>
void OutputBandPlot(Kpoint<KpointType> ** Kptr)
{

    double kx, ky, kz;
    char filename[4*MAX_PATH];
    FILE *bs_f;
    double *eig_all;
    int nspin = ct.spin_flag +1;

    int tot_num_eigs = nspin * ct.num_kpts * ct.num_states;
    eig_all = new double[tot_num_eigs];


    for(int idx = 0; idx < tot_num_eigs; idx++) eig_all[idx] = 0.0;

    for (int is = 0; is < ct.num_states; is++)
    {
        for(int ik = 0; ik < ct.num_kpts_pe; ik++)
        {
            int kpt = pct.kstart + ik;
            eig_all[pct.spinpe * ct.num_kpts * ct.num_states + is * ct.num_kpts + kpt] 
                = Kptr[ik]->Kstates[is].eig[0] * Ha_eV;
        }
    }

    GlobalSums (eig_all, tot_num_eigs, pct.kpsub_comm);


    double *x = new double[ct.num_kpts];
    x[0] = 0.0;
    for(int ik = 1; ik < ct.num_kpts; ik++)
    {
        kx = (ct.kp[ik].kvec[0]-ct.kp[ik-1].kvec[0]);
        ky = (ct.kp[ik].kvec[1]-ct.kp[ik-1].kvec[1]);
        kz = (ct.kp[ik].kvec[2]-ct.kp[ik-1].kvec[2]);
        x[ik] = x[ik-1] + sqrt(kx*kx + ky*ky + kz*kz);
    }

    if(pct.imgpe == 0)
    {

        for(int ispin = 0; ispin < nspin; ispin++)
        {
            snprintf(filename, sizeof(filename) - 1, "%s_spin%d%s", ct.basename, ispin, ".bandstructure.dat");
            bs_f = fopen (filename, "w");
            if(!bs_f) {
                rmg_printf("Unable to write band plot data.\n");
                return;
            }

            double eig_max = -1000.0;
            double eig_min = 1000.0;
            for (int is = 0; is < ct.num_states; is++)
            {
                for(int ik = 0; ik < ct.num_kpts; ik++)
                {
                    int idx = ispin * ct.num_kpts * ct.num_states + is * ct.num_kpts + ik;

                    fprintf (bs_f, "\n %f  %16.8f ", x[ik], eig_all[idx]);
                    eig_max = std::max(eig_max, eig_all[idx]);
                    eig_min = std::min(eig_min, eig_all[idx]);


                }

                fprintf (bs_f, "\n &&");
            }

            fclose (bs_f);


            eig_min -= 1.0;

            snprintf(filename, sizeof(filename) - 1, "%s_spin%d%s", ct.basename, ispin, ".bandstructure.xmgr");
            bs_f = fopen (filename, "w");
            if(!bs_f) {
                rmg_printf("Unable to write band plot data.\n");
                return;
            }
            fprintf (bs_f, "@version 50123\n");
            fprintf (bs_f, "@g0 on\n");
            fprintf (bs_f, "@with g0\n");
            fprintf (bs_f, "@    world  %f, %f, %f, %f\n", x[0], eig_min, x[ct.num_kpts-1], eig_max);
            fprintf (bs_f, "@    xaxis  on\n");
            fprintf (bs_f, "@    xaxis  tick spec type both\n");

            int num_tick = 0;
            for(int ik = 0; ik < ct.num_kpts; ik++)
            {
                if(strlen(ct.kp[ik].symbol))
                {
                    fprintf (bs_f, "@    xaxis  tick major %d, %f\n", num_tick, x[ik]);
                    fprintf (bs_f, "@    xaxis  ticklabel %d, \"%s\"\n", num_tick, ct.kp[ik].symbol);
                    num_tick++;
                }
            }
            fprintf (bs_f, "@    xaxis  tick spec %d\n", num_tick);
            fprintf (bs_f, "@    xaxis  ticklabel char size 1.75\n");

            fprintf(bs_f, "@    yaxis  on\n");
            fprintf(bs_f, "@    yaxis  label \"E(eV)\"\n");
            fprintf(bs_f, "@    yaxis  label char size 2.000000\n");
            fprintf(bs_f, "@    yaxis  tick on\n");
            fprintf(bs_f, "@    yaxis  tick major 5\n");
            fprintf(bs_f, "@    yaxis  tick minor ticks 1\n");
            fprintf(bs_f, "@    yaxis  tick major linewidth 2.0\n");
            fprintf(bs_f, "@    yaxis  tick minor linewidth 1.5\n");
            fprintf(bs_f, "@    yaxis  ticklabel char size 1.750000\n");
            fprintf(bs_f, "@    frame type 0\n");
            fprintf(bs_f, "@    frame linestyle 1\n");
            fprintf(bs_f, "@    frame linewidth 2.0\n");
            fprintf(bs_f, "@    frame color 1\n");
            fprintf(bs_f, "@    frame background color 0\n");

            for(int ik = 1; ik < ct.num_kpts-1; ik++)
            {
                if(strlen(ct.kp[ik].symbol))
                {
                    fprintf(bs_f, "@with line\n");
                    fprintf(bs_f, "@    line on\n");
                    fprintf(bs_f, "@    line loctype world\n");
                    fprintf(bs_f, "@    line g0\n");
                    fprintf(bs_f, "@    line %f, %f, %f, %f\n", x[ik], eig_min, x[ik], eig_max);
                    fprintf(bs_f, "@    line linewidth 1.0\n");
                    fprintf(bs_f, "@    line linestyle 1\n");
                    fprintf(bs_f, "@    line color 1\n");
                    fprintf(bs_f, "@line def\n");
                }
            }

            for (int is = 0; is < ct.num_states; is++)
            {
                fprintf(bs_f, "@    s%d hidden false\n", is);
                fprintf(bs_f, "@    s%d type xy\n", is);
                fprintf(bs_f, "@    s%d line type 1\n", is);
                fprintf(bs_f, "@    s%d line linestyle 1\n", is);
                fprintf(bs_f, "@    s%d line linewidth 2.0\n", is);
                fprintf(bs_f, "@    s%d line color 1\n", is);
            }


            for (int is = 0; is < ct.num_states; is++)
            {
                fprintf(bs_f, "@target G0.S%d\n", is);
                fprintf(bs_f, "@type xy\n");

                for(int ik = 0; ik < ct.num_kpts; ik++)
                {
                    int idx = ispin * ct.num_kpts * ct.num_states + is * ct.num_kpts + ik;

                    fprintf (bs_f, "\n %f  %16.8f ", x[ik], eig_all[idx]);

                }

                fprintf (bs_f, "\n &&");
            }

            fclose (bs_f);

        }


    }

}



