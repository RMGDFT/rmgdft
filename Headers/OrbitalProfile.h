#ifndef RMG_ORBITAL_PROFILE_H
#define RMG_ORBITAL_PROFILE_H 1

#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "transition.h"
#include <boost/algorithm/string.hpp>

extern std::unordered_map<std::string, std::string> SymbolToConfig;

// OrbitalProfile contains aggregate information on atomic orbitals
class OrbitalProfile
{

public:
    OrbitalProfile(void)
    {
        double acc = 0;
        for (int ion = 0; ion < ct.num_ions; ion++)
        {   
            /* Generate ion pointer */
            ION *iptr = &Atoms[ion];
            
            /* Get species type and electronic configuration */
            SPECIES *sp = &Species[iptr->species];
            std::string symbol(sp->atomic_symbol);
            std::string config = SymbolToConfig[symbol]; 
            std::vector<std::string> shells;
            boost::split(shells, config, boost::is_any_of(" "));

            int jx = 0;
            // Iterate from outer shells in.
            for (auto it = shells.rbegin(); it != shells.rend(); ++it)
            {
                std::string pp = *it;
                double shell_occ = std::stod(boost::algorithm::erase_head_copy(pp, 2));
                //if(pct.gridpe == 0) printf("OCC = %f\n", shell_occ);
                //if(pct.gridpe==0)std::cout << pp << std::endl;

                if (pp.find('s') != std::string::npos)
                {
                    if(shell_occ == 2.0) // fullshell
                    {
                        acc += 0.125;
                    }
                    else
                    {
                        double empty = 2 - shell_occ;
                        acc += 0.25 * empty;
                    }
                }
                else if(pp.find('p') != std::string::npos)
                {
                    if(shell_occ == 6.0)
                    {
                        acc += 0.125;
                    }
                    else
                    {
                        double empty = 6 - shell_occ;
                        acc += 0.125 * empty;
                    }
                }
                else if(pp.find('d') != std::string::npos)
                {
                    if(shell_occ == 10.0)
                    {
                        // Full d shell case then all symmetries exist
                        // in the valence band so just keep 20%.
                        acc += 1.0;
                    }
                    else
                    {
                        double empty = 10 - shell_occ;
                        // Partially occupied need to ensure all symmetries
                        // included. May be overkill for systems with a gap.
                        acc += 0.5 * empty;
                    }
                }
                else if(pp.find('f') != std::string::npos)
                {
                    if(shell_occ == 14.0)
                    {
                        acc += 1.4;
                    }
                    else
                    {
                        double empty = 14 - shell_occ;
                        acc += 0.7 * empty;
                    }
                }
                else
                {
                    throw RmgFatalException() <<  "Do you really have g-orbitals? " << 
                    __FILE__ << " at line " << __LINE__ << "\n";
                }
                jx++;
                if(jx == sp->num_atomic_waves) break;
            }
        }
        
        rec_unoccupied = std::ceil(acc);
        rec_unoccupied = std::max(5, rec_unoccupied);
        if(pct.gridpe==0 && ct.verbose)printf("rec_unoccupied = %d\n", rec_unoccupied);
   }
 
    int rec_unoccupied = 0;

};
#endif
