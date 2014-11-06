#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>



#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "MapElements.h"
#include "transition.h"


void ReadPseudo(int nspecies, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    // Default to norm conserving pseudopotentials. If any of the atomic species
    // are represented by ultrasoft then all of them must be handled as ultrasoft
    lc.norm_conserving_pp = true;

    for(int isp = 0;isp < nspecies;isp++) {
        SPECIES *sp = &lc.sp[isp];
        LoadUpf(sp);
        if(!sp->is_norm_conserving) lc.norm_conserving_pp = false;
        sp++; 
    } 


/**********************************************************************

   The UPF format includes fields to identify the XC type. Additionally
   the RMG input file has an option to specify the XC type to be used. The
   following rules are then used to determine the type of XC to be used
   in the calculation.

   If XC type is specified in the RMG input file it overides the value
   contained in the pseudopotentials.

   If there is no XC type specified in the RMG input file and all pseudopotentials
   have the same type then that type is used.

   If there is no XC type specified in the RMG input file and pseudopotentials
   specify more than one type of XC then throw an error and terminate since
   we have no way of knowing what the user intends.

**********************************************************************/

   // User specified the type explicity so we are done
   if(lc.xctype != AUTO_XC) return;
   
   std::string delims = " \t\n";

   bool uniform = true;
   std::set<std::string> short_names;

   for(int isp = 0;isp < nspecies;isp++) {

        SPECIES *sp = &lc.sp[isp];
        std::vector<std::string> func_components;
        std::string funcstr(sp->functional);
        boost::trim(funcstr);
        boost::algorithm::split( func_components, funcstr, boost::is_any_of(delims), boost::token_compress_on );

        // If there is only one entry then it must be a short name so insert it
        if(func_components.size() == 1) {

            std::string abbr = func_components.at(0);
            boost::to_upper(abbr); 
            short_names.emplace(abbr);

        }
        else {


        }
    

   }

}
