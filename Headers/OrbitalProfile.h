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
                }
                else if(ct.noncoll)
                {
                    totals[l] += 2*(2 * l + 1);
                }
                else
                {
                    totals[l] += 2 * l + 1;
                    occupied[l] += sp->atomic_wave_oc[ip] / 2;
                }
            }
        }

        for(size_t i=0;i < totals.size();i++)
            rec_unoccupied += ct.nspin * (totals[i] - occupied[i]) / 4;

        if(pct.gridpe == 0)
        //if(pct.gridpe == 0 && ct.verbose)
        {
            printf("Orbital Profile recomended unoccupied = %d\n", rec_unoccupied);
            printf("  S = %d  %d\n", totals[0], occupied[0]);
            printf("  P = %d  %d\n", totals[1], occupied[1]);
            printf("  D = %d  %d\n", totals[2], occupied[2]);
            printf("  F = %d  %d\n", totals[3], occupied[3]);
            printf("  G = %d  %d\n", totals[4], occupied[4]);
        }
   }
 
    std::vector<int> totals = {0, 0, 0, 0, 0};
    std::vector<int> occupied = {0, 0, 0, 0, 0};
    int rec_unoccupied = 0;

};
#endif
