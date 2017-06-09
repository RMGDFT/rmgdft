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

void ReadVelocities(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::string VelocityArray;
    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::vector<std::string> Velocities;
    std::unordered_map<std::string, InputKey *> NewMap;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("ionic_velocities", &VelocityArray, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "Ionic velocitys.\n",
                     "");
    
    
    If.LoadInputKeys();

    
    // Process velocities
    boost::trim(VelocityArray);
    boost::trim_if(VelocityArray, boost::algorithm::is_any_of("\"^"));

    if(VelocityArray.size() == 0) return;

    boost::algorithm::split( Velocities, VelocityArray, boost::is_any_of(line_delims), boost::token_compress_on );

    if(Velocities.size() != (size_t)lc.num_ions) {

        throw RmgFatalException() << "Ionic velocities must be present for " << lc.num_ions << " ions but only " << Velocities.size() << " found! Terminating.\n";

    }
    

    std::vector<std::string>::iterator it, it1;
    int nions = 0;
    for (it = Velocities.begin(); it != Velocities.end(); ++it) {

        std::string Velocity = *it;
        boost::trim(Velocity);

        std::vector<std::string> VelocityComponents;
        boost::algorithm::split( VelocityComponents, Velocity, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = VelocityComponents.size();
        if(ncomp != 3) throw RmgFatalException() << "Synax error in restart file near " << Velocity << "\n";

        it1 = VelocityComponents.begin();

        std::string xvelocity = *it1;
        boost::trim(xvelocity);
        lc.ions[nions].velocity[0] = std::atof(xvelocity.c_str());

        it1++;
        std::string yvelocity = *it1;
        boost::trim(yvelocity);
        lc.ions[nions].velocity[1] = std::atof(yvelocity.c_str());

        it1++;
        std::string zvelocity = *it1;
        boost::trim(zvelocity);
        lc.ions[nions].velocity[2] = std::atof(zvelocity.c_str());
        nions++;

    }

    if(nions > lc.num_ions)
        throw RmgFatalException() << "Inconsistency in number of ions: " << lc.num_ions << " was specified initially but " << nions << " were found.\n";

}
