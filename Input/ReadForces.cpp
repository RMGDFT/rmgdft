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



namespace Ri = RmgInput;

void ReadForces(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::string ForceArray;
    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::vector<std::string> Forces;
    std::unordered_map<std::string, InputKey *> NewMap;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("ionic_forces", &ForceArray, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "Ionic forces.\n",
                     "");
    
    
    If.LoadInputKeys();

    
    // Process forces
    boost::trim(ForceArray);
    boost::trim_if(ForceArray, boost::algorithm::is_any_of("\"^"));

    if(ForceArray.size() == 0) return;

    boost::algorithm::split( Forces, ForceArray, boost::is_any_of(line_delims), boost::token_compress_on );

    if(Forces.size() != (size_t)(4 * lc.num_ions)) {

        throw RmgFatalException() << "Ionic forces must be present for " << lc.num_ions << " ions but only " << Forces.size() << " found! Terminating.\n";

    }
    

    std::vector<std::string>::iterator it, it1;
    int nions = 0;
    int fpt = 0;   // We save the last 4 sets of forces in the restart file so have to read back that many too
    for (it = Forces.begin(); it != Forces.end(); ++it) {

        std::string Force = *it;
        boost::trim(Force);

        std::vector<std::string> ForceComponents;
        boost::algorithm::split( ForceComponents, Force, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = ForceComponents.size();
        if(ncomp != 3) throw RmgFatalException() << "Synax error in restart file near " << Force << "\n";

        it1 = ForceComponents.begin();

        std::string xforce = *it1;
        boost::trim(xforce);
        lc.ions[nions].force[fpt][0] = std::atof(xforce.c_str());

        it1++;
        std::string yforce = *it1;
        boost::trim(yforce);
        lc.ions[nions].force[fpt][1] = std::atof(yforce.c_str());

        it1++;
        std::string zforce = *it1;
        boost::trim(zforce);
        lc.ions[nions].force[fpt][2] = std::atof(zforce.c_str());
        nions++;
        // Check if we have read forces for all ions. If yes reset nions
        if(nions == lc.num_ions) {
            nions = 0;
            fpt++;
        }

    }

    if(nions > lc.num_ions)
        throw RmgFatalException() << "Inconsistency in number of ions: " << lc.num_ions << " was specified initially but " << nions << " were found.\n";

}
