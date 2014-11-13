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
#include "typedefs.h"



/**********************************************************************

    Reads internal RMG atomic format which consists of a single atom
    per line with the entire block of lines enclosed in double quotes
    and the block identified by an atoms keyword.

    Column 1: num of orbitals center at the following cnter
    Column 2: x
    Column 3: y
    Column 4: z  orbita center
    Column 5: radius
    Column 6: 
    Column 7: 
            

Example:

kpoints = "
    0.0  0.0  0.0   0.5
    0.5  0.5  0.5   0.5
"


    cfile        = name of the file containing the kpoint info
    InputMap     = Control Map. May not be needed by all atomic input
                   drivers but is useful for reading the RMG format.
    
**********************************************************************/

namespace Ri = RmgInput;

void ReadOrbitals(char *cfile, STATE  *states, int *state_to_ion, MPI_Comm comm)
{

    std::string OrbitalArray;
    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::vector<std::string> Orbitals;
    std::unordered_map<std::string, InputKey *> NewMap;
    int natom;
    int i, ni;
    double x, y, z,radius;
    double bohr;
    int num_tem, movable, frozen;


    RmgInputFile If(cfile, NewMap, comm);

    If.RegisterInputKey("orbitals", &OrbitalArray, "",
            CHECK_AND_FIX, REQUIRED,
            "orbitals list and their radius \n",
            "");


    If.LoadInputKeys();

    bohr = 1.0;
    if (Verify ("length_units", "Angstrom", NewMap))
    {
        bohr= 1.88972612499359;
    }

    // Process atoms
    boost::trim(OrbitalArray);
    boost::trim_if(OrbitalArray, boost::algorithm::is_any_of("\"^"));

    boost::algorithm::split( Orbitals, OrbitalArray, boost::is_any_of(line_delims), boost::token_compress_on );


    std::vector<std::string>::iterator it, it1;
    natom=0;
    ni = 0;

    for (it = Orbitals.begin(); it != Orbitals.end(); ++it) {

        std::string Orbital = *it;
        boost::trim(Orbital);

        std::vector<std::string> OrbitalComponents;
        boost::algorithm::split( OrbitalComponents, Orbital, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = OrbitalComponents.size();
        if((ncomp < 7) ) throw RmgFatalException() << "Synax error in orbital information near " << Orbital << "\n";

        it1 = OrbitalComponents.begin();

        std::string nstr = *it1;
        num_tem = std::atoi(nstr.c_str());


        it1++;
        std::string xstr = *it1;
        x = std::atof(xstr.c_str()) * bohr;
        it1++;
        std::string ystr = *it1;
        y = std::atof(ystr.c_str()) * bohr;
        it1++;
        std::string zstr = *it1;
        z = std::atof(zstr.c_str()) * bohr;
        it1++;
        std::string rstr = *it1;
        radius = std::atof(rstr.c_str());

        it1++;
        rstr = *it1;
        movable = std::atoi(rstr.c_str());
        it1++;
        rstr = *it1;
        frozen = std::atoi(rstr.c_str());


        for(i = 0; i < num_tem; i++)
        {
            states[ni + i].crds[0] = x * bohr;
            states[ni + i].crds[1] = y * bohr;
            states[ni + i].crds[2] = z * bohr;
            states[ni + i].radius = radius;
            states[ni + i].movable = movable;
            states[ni + i].frozen = frozen;
            states[ni + i].n_orbital_same_center = num_tem;
            states[ni + i].gaussian_orbital_index = i;
            state_to_ion[ni+i] = natom;
            states[ni+i].atomic_orbital_index = i;
        }


        natom++;
        ni += num_tem;

    }


}
