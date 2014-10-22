#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;
#include <iostream> 
#include <fstream>
#include <sstream>
#include <iterator>
#include <string> 
#include <cfloat> 
#include <climits> 
#include <unordered_map>
#include <set>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include "MapElements.h"
#include "transition.h"
#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "CheckValue.h"
#include "RmgException.h"
#include "RmgInputFile.h"
#include "InputOpts.h"


/**********************************************************************

    For detailed documentation on how the input system works look at
    Input/ReadCommon.cpp. This function is used to read specific
    information related to ionic positions, molecular dynamics and relaxations.

**********************************************************************/

namespace Ri = RmgInput;

void ReadDynamics(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::set<std::string> SpeciesTypes;
    std::list<std::string> IonSpecies;

    // Determine input file format and dispatch to the correct driver routine
    ReadRmgAtoms(cfile, SpeciesTypes, IonSpecies, lc, InputMap);

    // Atoms have been read, now take care of conversion to internal coordinates
    if (Verify ("atomic_coordinate_type", "Cell Relative", InputMap)) {

        for(int ion = 0;ion < lc.num_ions;ion++) {

            lc.ions[ion].xtal[0] = lc.ions[ion].crds[0];
            lc.ions[ion].xtal[1] = lc.ions[ion].crds[1];
            lc.ions[ion].xtal[2] = lc.ions[ion].crds[2];
            to_cartesian(lc.ions[ion].xtal, lc.ions[ion].crds);

        }

    }
    else if(Verify ("crds_units", "Angstrom", InputMap)) {

        for(int ion = 0;ion < lc.num_ions;ion++) {

            lc.ions[ion].crds[0] *= A_a0;
            lc.ions[ion].crds[1] *= A_a0;
            lc.ions[ion].crds[2] *= A_a0;

        }

    }

    // SpeciesType holds the number of species found
    lc.num_species = SpeciesTypes.size();
    ct.sp = new SPECIES[ct.num_species];

    int isp = 0;
    for(auto it = SpeciesTypes.begin();it != SpeciesTypes.end(); ++it) {
        std::string AtomicSymbol = *it;
        lc.sp[isp].atomic_symbol = new char[4];
        std::strncpy(lc.sp[isp].atomic_symbol, AtomicSymbol.c_str(), 4);
        std::strncpy(lc.sp[isp].pseudo_symbol, AtomicSymbol.c_str(), 4);
        std::strncpy(lc.sp[isp].pseudo_filename, "./@Internal", sizeof(lc.sp[isp].pseudo_filename));
        LoadUpf(&ct.sp[isp]);
        isp++;
    }

    // Assign species type for each ion
    int species = 0;
    for(auto it = SpeciesTypes.begin();it != SpeciesTypes.end(); ++it) {

        std::string AtomicSymbol1 = *it;
        int ion = 0;
        for(auto it1 = IonSpecies.begin();it1 != IonSpecies.end(); ++it1) {

            std::string AtomicSymbol2 = *it1;
            if(!AtomicSymbol1.compare(AtomicSymbol2)) {
                lc.ions[ion].species = species;
            }
            ion++;
            
        }
        species++;

    }

}
