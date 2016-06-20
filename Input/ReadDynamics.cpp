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

    std::string string_tem;

    std::unordered_map<std::string, InputKey *> NewMap;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("kpoint_distribution", &pct.pe_kpoint, -INT_MAX, INT_MAX, -1,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "");

    If.LoadInputKeys();


    // Read atoms
    ReadRmgAtoms(cfile, SpeciesTypes, IonSpecies, lc, InputMap);

    ReadTFAtoms(cfile, SpeciesTypes, IonSpecies, lc, InputMap);

    // Forces and velocities (if present)
    ReadForces(cfile, lc, InputMap);
    ReadVelocities(cfile, lc, InputMap);

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
        lc.sp[isp].atomic_symbol = new char[4]();
        std::strncpy(lc.sp[isp].atomic_symbol, AtomicSymbol.c_str(), 3);
        std::strncpy(lc.sp[isp].pseudo_filename, "./@Internal", sizeof(lc.sp[isp].pseudo_filename));
        isp++;
    }

    // Next determine if an external pseudopotential has been specified that overrides the internal one.
    InputKey *Key;
    try {
        std::vector<std::string> lines;
        std::string delims = "\r\n^";
        std::string field_delims = " \t";
        Key = InputMap.at("pseudopotential");
        boost::algorithm::split( lines, Key->Readstr, boost::is_any_of(delims), boost::token_compress_on );
        std::vector<std::string>::iterator it;
        for (it = lines.begin(); it != lines.end(); ++it) {
            std::string pline = *it;
            boost::trim_if(pline, boost::algorithm::is_any_of("\" \t"));
            std::vector<std::string> fields;
            boost::algorithm::split( fields, pline, boost::is_any_of(field_delims), boost::token_compress_on );
            if(fields.size() == 2) {

                // Search the species structure for a matching symbol
                boost::trim_if(fields[0], boost::algorithm::is_any_of("\" \t"));
                for(int isp = 0;isp < lc.num_species;isp++) {
                    if(!std::strcmp(fields[0].c_str(), lc.sp[isp].atomic_symbol)) {

                        string_tem = std::string(pct.image_path[pct.thisimg]) + fields[1];  
                        std::strncpy(lc.sp[isp].pseudo_filename, string_tem.c_str(), sizeof(lc.sp[isp].pseudo_filename));

                    }
                }

            }

        }

    }
    catch (const std::out_of_range& oor) {
        // no pseudpopotential tag which means only use internals
    }

    // Load pseudopotentials
    for(int isp = 0;isp < lc.num_species;isp++) {
        LoadUpf(&ct.sp[isp]);
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

