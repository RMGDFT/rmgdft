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
#include <boost/filesystem.hpp>
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


#if OPENBABEL_LIBS
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#endif


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
    IonSpecies   = a std::list that contains the atomic symbol of each
                   ion in the order in which they were read from cfile.
    InputMap     = Control Map. May not be needed by all atomic input
                   drivers but is useful for reading the RMG format.
    
**********************************************************************/

namespace Ri = RmgInput;

void ReadTFAtoms(char *cfile, std::set<std::string>& SpeciesTypes, std::list<std::string>& IonSpecies, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::string AtomArray;
    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::vector<std::string> Atoms;
    std::unordered_map<std::string, InputKey *> NewMap;
    int nions = 0;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("tf_atoms", &AtomArray, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "Positions and gaussian charge parameters for simplified waters .\n",
                     "");

    
    
    If.LoadInputKeys();


    // Process atoms
    boost::trim(AtomArray);
    boost::trim_if(AtomArray, boost::algorithm::is_any_of("\"^"));

    if(AtomArray.size() == 0) 
    {
	lc.num_tfions = 0;
	//throw RmgFatalException() << "1Found" << lc.num_tfions << " TF ions \n";
	return;
    }


    boost::algorithm::split( Atoms, AtomArray, boost::is_any_of(line_delims), boost::token_compress_on );



    lc.num_tfions = Atoms.size();
    printf("\n Found %d TF ions", lc.num_tfions);
    ///throw RmgFatalException() << "Found" << lc.num_tfions << " TF ions \n";
    lc.tf_ions = new TF_ION[lc.num_tfions]();

    std::vector<std::string>::iterator it, it1;
    for (it = Atoms.begin(); it != Atoms.end(); ++it) {

        if(nions > lc.num_tfions)
            throw RmgFatalException() << "Inconsistency in number of TF ions: " << lc.num_tfions << " was specified initially but " << nions << " were found.\n";

        
	std::string Atom = *it;
        boost::trim(Atom);

        std::vector<std::string> AtomComponents;
        boost::algorithm::split( AtomComponents, Atom, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = AtomComponents.size();
        if((ncomp < 8) || (ncomp > 8)) throw RmgFatalException() << "Syntax error in tf_ions, 8 arguments are needed on each line, but" << ncomp << "were found \n";


        
	// First field should be an atomic symbol
        it1 = AtomComponents.begin();
        std::string sp = *it1;
        boost::trim(sp);
        
	// Valid atomic symbol? GetAtomicMass will throw a fatal exception if the symbol is not valid.
        GetAtomicMass(sp);
        
        
	//Process coordinates and other data on the line
        it1++;
	std::string xstr = *it1;
        boost::trim(xstr);
	lc.tf_ions[nions].crds[0] = std::atof(xstr.c_str());


	it1++;
	std::string ystr = *it1;
        boost::trim(ystr);
	lc.tf_ions[nions].crds[1] = std::atof(ystr.c_str());
	
	it1++;
	std::string zstr = *it1;
        boost::trim(zstr);
	lc.tf_ions[nions].crds[2] = std::atof(zstr.c_str());
	
	it1++;
	std::string qstr = *it1;
        boost::trim(qstr);
	lc.tf_ions[nions].q = std::atof(qstr.c_str());
	
	it1++;
	std::string alphastr = *it1;
        boost::trim(alphastr);
	lc.tf_ions[nions].alpha = std::atof(alphastr.c_str());
	
	it1++;
	std::string q0str = *it1;
        boost::trim(q0str);
	lc.tf_ions[nions].q0 = std::atof(q0str.c_str());
	
	it1++;
	std::string alpha0str = *it1;
        boost::trim(alpha0str);
	lc.tf_ions[nions].alpha0 = std::atof(alpha0str.c_str());
      


        nions++;

    }

    if(nions > lc.num_tfions)
        throw RmgFatalException() << "Inconsistency in number of ions: " << lc.num_tfions << " was specified initially but " << nions << " were found.\n";

    //printf("\n tf_atoms read OK. Number is %d", lc.num_tfions);

}
