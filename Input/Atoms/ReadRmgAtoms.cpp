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
#include <list>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "MapElements.h"
#include "transition.h"
#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgInputFile.h"
#include "InputOpts.h"


/**********************************************************************

    Reads internal RMG atomic format which consists of a single atom
    per line with the entire block of lines enclosed in double quotes
    and the block identified by an atoms keyword.

    Column 1: Atomic symbol
    Column 2: x-coordinate
    Column 3: y-coordinate
    Column 4: z-coordinate
    Column 5: boolean value indicating whether the ion is movable or not.
              can be omitted and the default is movable(1).

Example:

atoms = "
N    15.3293    12.9600    19.2070   1
C    17.5363    12.9600    17.7857    1
H    13.7313    15.1732    19.2251    1
O    14.4083    17.2893    17.8251    1
C    11.1292    14.3262    19.2527    1
C     9.3055    15.6378    17.8881    1
"

    The coordinates may be either cell relative or absolute in units
    of either bohr or Angstroms. The driver routine is not responsible
    for transforming the input coordinates into the internal format
    which is handled at a higher level.
    
    Individual driver routines for different formats require 3 arguments.

    cfile        = name of the file containing the ionic information.
    SpeciesTypes = a std::set that is empty on entry and contains the
                   atomic symbol of each unique species on exit.
    Species      = a std::list that contains the atomic symbol of each
                   ion in the order in which they were read from cfile.
    InputMap     = Control Map. May not be needed by all atomic input
                   drivers but is useful for reading the RMG format.
    
**********************************************************************/

namespace Ri = RmgInput;

void ReadRmgAtoms(char *cfile, std::set<std::string>& SpeciesTypes, std::list<std::string>& Species, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::string AtomArray;
    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::vector<std::string> Atoms;
    int nions = 0;

    RmgInputFile If(cfile, InputMap);

    If.RegisterInputKey("atoms", &AtomArray, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "Ionic species and coordinate information.\n",
                     "");
    
    
    If.LoadInputKeys();

    // Process atoms
    boost::trim(AtomArray);
    boost::trim_if(AtomArray, boost::algorithm::is_any_of("\"^"));

    boost::algorithm::split( Atoms, AtomArray, boost::is_any_of(line_delims), boost::token_compress_on );
    lc.num_ions = Atoms.size();
    lc.ions = new ION[lc.num_ions]();

    std::vector<std::string>::iterator it, it1;
    for (it = Atoms.begin(); it != Atoms.end(); ++it) {

        if(nions > lc.num_ions)
            throw RmgFatalException() << "Inconsistency in number of ions: " << lc.num_ions << " was specified initially but " << nions << " were found.\n";

        std::string Atom = *it;
        boost::trim(Atom);

        std::vector<std::string> AtomComponents;
        boost::algorithm::split( AtomComponents, Atom, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = AtomComponents.size();
        if((ncomp < 4) || (ncomp > 5)) throw RmgFatalException() << "Synax error in ionic information near " << Atom << "\n";

        // First field should be an atomic symbol
        it1 = AtomComponents.begin();
        std::string sp = *it1;
        boost::trim(sp);

        // Valid atomic symbol? GetAtomicMass will throw a fatal exception if the symbol is not valid.
        GetAtomicMass(sp);

        // Is valid so make an entry in the SpeciesTypes set and in the Species list. SpeciesTypes set only contains
        // one entry for each unique species type while Species contains lc.num_ions entries (one for each ion)
        SpeciesTypes.emplace(sp);
        Species.emplace_back(sp);
      
        // Look for the coordinates
        it1++;
        std::string xstr = *it1;
        lc.ions[nions].crds[0] = std::atof(xstr.c_str());
        it1++;
        std::string ystr = *it1;
        lc.ions[nions].crds[1] = std::atof(ystr.c_str());
        it1++;
        std::string zstr = *it1;
        lc.ions[nions].crds[2] = std::atof(zstr.c_str());

        int movable = 1;
        std::string smov;
        if(ncomp == 5) {
            it1++;
            smov = *it1;
            boost::trim(smov);
            if( !smov.compare("0") ) movable = 0;
            if( !smov.compare("no") ) movable = 0;
            if( !smov.compare("No") ) movable = 0;
            if( !smov.compare("NO") ) movable = 0;
        }

        lc.ions[nions].movable = movable;

        nions++;

    }

}
