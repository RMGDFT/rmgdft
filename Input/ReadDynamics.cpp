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

typedef struct {double U; int index; std::string label; double J[3]; } HUBBARD_INFO;
void parse_Hubbard_info(std::string what, std::vector<HUBBARD_INFO> &hinfo, std::unordered_map<std::string, InputKey *>& InputMap, CONTROL& lc);

namespace Ri = RmgInput;

void ReadDynamics(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::set<std::string> SpeciesTypes;
    std::list<std::string> IonSpecies;

    std::unordered_map<std::string, InputKey *> NewMap;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("kpoint_distribution", &pct.pe_kpoint, -INT_MAX, INT_MAX, -1,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "");
    If.RegisterInputKey("FermiEnergy", &lc.efermi, -1.0e9, 1.0e9, 0.0,
            CHECK_AND_FIX, OPTIONAL,
            "read in Fermi energy for band strcuture plot ",
            "Only needed for band structure  ", MISC_OPTIONS);


    If.LoadInputKeys();


    // Read atoms
    ReadRmgAtoms(cfile, SpeciesTypes, IonSpecies, Atoms, lc);

    ReadTFAtoms(cfile, SpeciesTypes, IonSpecies, lc, InputMap);

    // Forces and velocities (if present)
    ReadForces(cfile, lc, InputMap);
    ReadVelocities(cfile, lc, InputMap);

    // SpeciesType holds the number of species found
    lc.num_species = SpeciesTypes.size();
    Species.resize(lc.num_species);

    std::string nums = "-.0123456789";
    int isp = 0;
    for(auto it = SpeciesTypes.begin();it != SpeciesTypes.end(); ++it) {
        std::string AtomicSymbol = *it;
        for (char c: nums) AtomicSymbol.erase(std::remove(AtomicSymbol.begin(), AtomicSymbol.end(), c), AtomicSymbol.end());
        Species[isp].atomic_symbol = new char[4]();
        std::strncpy(Species[isp].atomic_symbol, AtomicSymbol.c_str(), 3);
        Species[isp].pseudo_filename = std::string("./@Internal");
        Species[isp].index = isp;
        isp++;
    }

    // Next determine if an external pseudopotential has been specified that overrides the internal one.
    try {
        InputKey *Key;
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
                    if(!std::strcmp(fields[0].c_str(), Species[isp].atomic_symbol))
                    {
                        Species[isp].pseudo_filename = std::string(pct.image_path[pct.thisimg]) + fields[1];  
                    }
                }
            }
        }
    }
    catch (const std::out_of_range& oor) {
        // no pseudpopotential tag which means only use internals
    }

    for(int isp = 0;isp < lc.num_species;isp++)
    {
        Species[isp].Hubbard_U = 0.0;
        Species[isp].Hubbard_J[0] = 0.0;
        Species[isp].Hubbard_J[1] = 0.0;
        Species[isp].Hubbard_J[2] = 0.0;
    }

    std::vector<HUBBARD_INFO> hinfo;
    parse_Hubbard_info("Hubbard_U", hinfo, InputMap, lc);
    for(auto it=hinfo.begin();it != hinfo.end();++it) 
    {
        HUBBARD_INFO h = *it;
        Species[h.index].Hubbard_U = h.U;
        Species[h.index].ldaU_label = h.label;
        Species[h.index].is_ldaU = true;
        Species[h.index].Hubbard_J[0] = h.J[0];
        Species[h.index].Hubbard_J[1] = h.J[1];
        Species[h.index].Hubbard_J[2] = h.J[2];
    
    }


    // Load pseudopotentials
    for(int isp = 0;isp < lc.num_species;isp++)
    {
        LoadPseudo(&Species[isp]);
    }

    // Assign species type for each ion
    int species = 0;
    for(auto it = SpeciesTypes.begin();it != SpeciesTypes.end(); ++it) {

        std::string AtomicSymbol1 = *it;
        int ion = 0;
        for(auto it1 = IonSpecies.begin();it1 != IonSpecies.end(); ++it1) {

            std::string AtomicSymbol2 = *it1;
            if(!AtomicSymbol1.compare(AtomicSymbol2)) {
                Atoms[ion].species = species;
                Atoms[ion].Type = &Species[species];
            }
            ion++;

        }
        species++;

    }

}


// Used to parse out Hubbard_U and Hubbard_J vals
void parse_Hubbard_info(std::string what, std::vector<HUBBARD_INFO> &hinfo, std::unordered_map<std::string, InputKey *>& InputMap, CONTROL& lc)
{
    // Extract Hubbard U info if present
    try {
        InputKey *Key;
        std::vector<std::string> lines;
        std::string delims = "\r\n^";
        std::string field_delims = " \t";
        Key = InputMap.at(what);
        boost::algorithm::split( lines, Key->Readstr, boost::is_any_of(delims), boost::token_compress_on );
        std::vector<std::string>::iterator it;
        for (it = lines.begin(); it != lines.end(); ++it) {
            std::string pline = *it;
            boost::trim_if(pline, boost::algorithm::is_any_of("\" \t"));
            std::vector<std::string> fields;
            boost::algorithm::split( fields, pline, boost::is_any_of(field_delims), boost::token_compress_on );
            if(fields.size() >= 2)
            {
                // Search the species structure for a matching symbol
                boost::trim_if(fields[0], boost::algorithm::is_any_of("\" \t"));
                for(int isp = 0;isp < lc.num_species;isp++)
                {
                    // Found so try to extract a parameter
                    if(!std::strcmp(fields[0].c_str(), Species[isp].atomic_symbol))
                    {
                        HUBBARD_INFO h;
                        h.U = std::stod(fields[1]) / Ha_eV;
                        h.index = isp;
                        //h.label = "3d"; //default d orbital
                        h.J[0] = 0.0;
                        h.J[1] = 0.0;
                        h.J[2] = 0.0;

                        if(fields.size() < 3)
                            throw RmgFatalException() << "Missing orbital label in Hubbard_U specification. The label must match the orbital label from the pseudopotential. Terminating.\n";

                        h.label = fields[2];
                        if(fields.size() >= 4) h.J[0] = std::stod(fields[3]) / Ha_eV;
                        if(fields.size() >= 5) h.J[1] = std::stod(fields[4]) / Ha_eV;
                        if(fields.size() >= 6) h.J[2] = std::stod(fields[5]) / Ha_eV;
                        hinfo.push_back(h);
                    }
                }
            }
        }
    }
    catch (const std::out_of_range& oor) {
        // no Hubbard_U tag which is OK
    }
}

