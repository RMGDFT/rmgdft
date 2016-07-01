#include "portability.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>  
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
#include "RmgException.h"
#include "MapElements.h"
#include "InputOpts.h"
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

   std::vector<std::string> short_names;

   for(int isp = 0;isp < nspecies;isp++) {

        SPECIES *sp = &lc.sp[isp];
        std::vector<std::string> func_components;
        std::string funcstr(sp->functional);
        boost::trim_if(funcstr, boost::algorithm::is_any_of(" \t"));
        boost::to_upper(funcstr);
        boost::algorithm::split( func_components, funcstr, boost::is_any_of(delims), boost::token_compress_on );

        // If there is only one entry then it must be a short name so insert it
        if(func_components.size() == 1) {

            std::string abbr = func_components.at(0);
            boost::to_upper(abbr); 
            short_names.push_back(abbr);

        }
        else {

            // Must be a long name so search for the ones we support by
            // combining components into a single uppper case string and looking for a match
            std::string nstr;
            for(auto it = func_components.begin();it != func_components.end(); ++it) {
                boost::trim_if(*it, boost::algorithm::is_any_of(" \t"));
                nstr = nstr + *it;
            }
            boost::to_upper(nstr); 
            if(!nstr.compare("SLAPZ")) {
                short_names.push_back("PZ");
            }
            else if(!nstr.compare("SLAPZNOGXNOGC")) {
                short_names.push_back("PZ");
            }
            else if(!nstr.compare("SLAPWPBXPBC")) {
                short_names.push_back("PBE");
            }
            else if(!nstr.compare("SLAPWPBEPBE")) {
                short_names.push_back("PBE");
            }
            else if(!nstr.compare("SLAB88LYPBLYP")) {
                short_names.push_back("BLYP");
            }
            else if(!nstr.compare("B88P86")) {
                short_names.push_back("BP");
            }
            else if(!nstr.compare("SLAPWGGXGGC")) {
                short_names.push_back("PW91");
            }
            else if(!nstr.compare("LDA")) {
                short_names.push_back("PZ");
            }
            else {
                throw RmgFatalException() << "Unknown XC type in " << __FILE__ << " at line " << __LINE__ << ". Terminating.\n";
            }
            //std::cout << "FFF " << nstr << "LLL " << short_names.size() << std::endl;
            
        }

        // Remove any duplicates
        std::sort(short_names.begin(), short_names.end());
        std::vector<std::string>::iterator it;
        it = std::unique(short_names.begin(), short_names.end());
        short_names.resize(std::distance(short_names.begin(),it) );

        if(short_names.size() > 1) {
            throw RmgFatalException() << "Error: the pseudopotentials specified in the calculation have different exchange correlation types.\nEither switch to pseudopotentials with identical exchange correlation types or set a manual override in your input file using the exchange_correlation_type option.\n";
        }
        if(short_names.size() == 0) {
            throw RmgFatalException() << "No supported exchange_correlation_type found in the specified pseudopotentials. Try setting a manual override in your input file if you wish to use these pseudopotentials.\n";
        }

        // Finally classify
        std::string short_name = short_names.at(0);
        if(!short_name.compare("PZ")) {
            lc.xctype = LDA_PZ81; 
        }
        else if(!short_name.compare("PBE")) {
            lc.xctype = GGA_PBE;
        }
        else if(!short_name.compare("PBE")) {
            lc.xctype = GGA_PBE;
        }
        else if(!short_name.compare("BLYP")) {
            lc.xctype = GGA_BLYP;
        }
        else if(!short_name.compare("BP")) {
            lc.xctype = GGA_XB_CP;
        }
        else if(!short_name.compare("PW91")) {
            lc.xctype = GGA_XP_CP;
        }
        
   }

}
