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


#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"







/* Writes out the positions of the ions and the current forces on them */
void write_force (void)
{
    int num_movable = 0, maxfx_ion = 0, maxfy_ion = 0, maxfz_ion = 0, maxf_ion = 0;
    int num_movable_x = 0, num_movable_y = 0, num_movable_z = 0; 
    double avfx = 0.0, avfy = 0.0, avfz = 0.0, maxfx = 0.0, maxfy = 0.0, maxfz = 0.0;
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    double avf = 0.0;
    double maxf = 0.0, max_all_f = 0.0;
    double f2;
    double *fp;
    double efactor = ct.energy_output_conversion[ct.energy_output_units];
    const char *eunits = ct.energy_output_string[ct.energy_output_units].c_str();

    printf ("\n\n\n  IONIC POSITIONS [a0] AND FORCES [%s/a0]", eunits);
    printf     ("\n  Charge analysis using: "); 
    
    switch (ct.charge_analysis_type)
    {
	case CHARGE_ANALYSIS_NONE:
	    printf("No charge analysis performed");
	    break;
	
	case CHARGE_ANALYSIS_VORONOI:
	    printf("Voronoi Deformation Density");
	    break;

	default :
	    printf("Invalid Charge Analysis\n" );
    }

    printf("\n\n");

//    if (verify ("atom_constraints", NULL))
//    {
//        printf ("  CONSTRAINTS ON FORCES HAVE BEEN IMPOSED:\n\n");
//    }
//

    /*If no charge analysis has been performed set partial charges to 0*/
    if (ct.charge_analysis_type == CHARGE_ANALYSIS_NONE)
    {
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
	    Atoms[ion].partial_charge = 0.0;
    }

    printf
        ("@ION  Ion  Species       X           Y           Z       Charge       FX          FY         FZ      Movable\n");

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        ION &Atom = Atoms[ion];
        SPECIES &AtomType = Species[Atom.species];

        fp = Atom.force[ct.fpt[0]];

        printf ("@ION  %3lu  %4s     %10.7f  %10.7f  %10.7f   %6.3f   %10.7f  %10.7f  %10.7f  %1d %1d %1d\n",
                ion + 1,
                AtomType.atomic_symbol,
                Atom.crds[0], Atom.crds[1], Atom.crds[2], 
                 Atom.partial_charge,
		 efactor*fp[0], efactor*fp[1], efactor*fp[2], Atom.movable[0], Atom.movable[1], Atom.movable[2]);

        f2 = 0.0;
        if (Atom.movable[0])
        {

            num_movable_x++;

            avfx += fabs (fp[0]);


            if (fabs (fp[0]) > maxfx)
            {
                maxfx = fabs (fp[0]);
                maxfx_ion = ion;
            }

            f2 += fp[0] * fp[0];

        }

        if (Atom.movable[1])
        {
            avfy += fabs (fp[1]);
            if (fabs (fp[1]) > maxfy)
            {
                maxfy = fabs(fp[1]);
                maxfy_ion = ion;
            }
            f2 += fp[1] * fp[1];

        }

        if (Atom.movable[2])
        {
            avfz += fabs (fp[2]);
            if (fabs (fp[2]) > maxfz)
            {
                maxfz = fabs(fp[2]);
                maxfz_ion = ion;
            }
            f2 += fp[2] * fp[2];

        }





        if (f2 > maxf)
        {
            maxf = f2;
            maxf_ion = ion;
        } 

        avf += f2;

        sumx += fp[0];
        sumy += fp[1];
        sumz += fp[2];
    }



    num_movable = num_movable_x + num_movable_y + num_movable_z;
    if (num_movable_x != 0)
        avfx /= num_movable_x;
    if (num_movable_y != 0)
        avfy /= num_movable_y;
    if (num_movable_z != 0)
        avfz /= num_movable_z;

    if (num_movable != 0)
    {
        maxf = sqrt (maxf);
        avf = sqrt (avf / num_movable);
        max_all_f = rmg_max (maxfx, maxfy);
        max_all_f = rmg_max (max_all_f, maxfz);


        printf ("\n");
        printf (" mean FX      = %12.8f %s/a0\n", efactor*avfx, eunits);
        printf (" mean FY      = %12.8f %s/a0\n", efactor*avfy, eunits);
        printf (" mean FZ      = %12.8f %s/a0\n", efactor*avfz, eunits);

        printf ("\n");
        printf (" max FX       = %12.8f %s/a0   (ion %d)\n", efactor*maxfx, eunits, maxfx_ion + 1);
        printf (" max FY       = %12.8f %s/a0   (ion %d)\n", efactor*maxfy, eunits, maxfy_ion + 1);
        printf (" max FZ       = %12.8f %s/a0   (ion %d)\n", efactor*maxfz, eunits, maxfz_ion + 1);
        printf (" max F[x,y,z] = %12.8f %s/a0\n", efactor*max_all_f, eunits);
        printf (" max |F|      = %12.8f %s/a0   (ion %d)\n", efactor*maxf, eunits, maxf_ion + 1);
        if (ct.forceflag == MD_FASTRLX)
        {
            progress_tag ();
            printf (" tolerance    = %12.8f %s/a0\n", efactor*ct.thr_frc, eunits);
        }

        printf ("\n");
        if(!ct.renormalize_forces)
        {

            printf (" sum FX       = %12.8f %s/a0\n", efactor*sumx, eunits);
            printf (" sum FY       = %12.8f %s/a0\n", efactor*sumy, eunits);
            printf (" sum FZ       = %12.8f %s/a0\n", efactor*sumz, eunits);
            printf (" Average      = %12.8f %s/a0\n", efactor*(fabs (sumx) + fabs (sumy) + fabs (sumz)) / 3.0, eunits);
            printf (" sqrt < F^2 > = %12.8f %s/a0\n", efactor*avf, eunits);
            printf ("\n");

        }

    }


}                               /* end write_force */

/******/
