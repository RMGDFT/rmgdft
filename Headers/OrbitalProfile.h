#ifndef RMG_ORBITAL_PROFILE_H
#define RMG_ORBITAL_PROFILE_H 1

#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "transition.h"

// OrbitalProfile contains aggregate information on atomic orbitals
class OrbitalProfile
{

public:
    OrbitalProfile(void)
    {
        for (int ion = 0; ion < ct.num_ions; ion++)
        {   
            /* Generate ion pointer */
            ION *iptr = &Atoms[ion];
            
            /* Get species type */
            SPECIES *sp = &Species[iptr->species];
            
            for (int ip = 0; ip < sp->num_atomic_waves; ip++)
            {
                int l = sp->atomic_wave_l[ip];
                double jj = sp->atomic_wave_j[ip];
                if( sp->is_spinorb )
                {
                   totals[l] += (int)(2 * jj + 1 + 0.1);
                   total_atomic += (int)(2 * jj + 1 + 0.1);
                }
                else if(ct.noncoll)
                {
                    totals[l] += 2*(2 * l + 1);
                    total_atomic += 2*(2 * l + 1);
                }
                else
                {
                    totals[l] += 2 * l + 1;
                    total_atomic += 2 * l + 1;
                    occupied[l] += sp->atomic_wave_oc[ip];
                    total_occupied += sp->atomic_wave_oc[ip] / 2.0;
                }
            }
        }
        
#if 1
        for(size_t l=0;l < totals.size();l++)
        {
            switch(l)
            {
               case 0:
                   rec_unoccupied += std::floor(0.125*((double)totals[l] - std::floor(occupied[l]/2.0)));
                   break;
               case 1:
                   rec_unoccupied += std::floor(0.25*((double)totals[l] - std::floor(occupied[l]/2.0)));
                   break;
               default:
                   rec_unoccupied += std::floor(0.5*((double)totals[l] - std::floor(occupied[l]/2.0)));
            }
        }
        rec_unoccupied *= ct.nspin;
#endif
        if(pct.gridpe == 0 && ct.verbose)
        {
            printf("Orbital Profile recomended unoccupied = %d  %d  %f\n",
                    total_atomic, rec_unoccupied, total_occupied);
            printf("  S = %d  %f\n", totals[0], occupied[0]);
            printf("  P = %d  %f\n", totals[1], occupied[1]);
            printf("  D = %d  %f\n", totals[2], occupied[2]);
            printf("  F = %d  %f\n", totals[3], occupied[3]);
            printf("  G = %d  %f\n", totals[4], occupied[4]);
        }
   }
 
    std::vector<int> totals = {0, 0, 0, 0, 0};
    std::vector<double> occupied = {0, 0, 0, 0, 0};
    int rec_unoccupied = 0;
    int total_atomic = 0;
    double total_occupied = 0;

};
#endif
