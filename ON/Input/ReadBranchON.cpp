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

#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgInputFile.h"
#include "InputOpts.h"



/**********************************************************************

    read other stuff for ON branch

    cfile        = name of the file containing the kpoint info
    InputMap     = Control Map. May not be needed by all atomic input
                   drivers but is useful for reading the RMG format.
    
**********************************************************************/

namespace Ri = RmgInput;

void ReadBranchON(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::unordered_map<std::string, InputKey *> NewMap;

    std::string RMGfile;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    


    If.RegisterInputKey("start_mode", NULL, &lc.runflag, "",
                     CHECK_AND_TERMINATE, REQUIRED, start_mode,
                     "Type of run. Choices are \"Random Start\", \"Restart From File\", \"FIREBALL Start\", \"Gaussian Start\", or \"Restart TDDFT\".\n", 
                     "start_mode must be one of  \"Random Start\", \"Restart From File\", \"FIREBALL Start\", \"Gaussian Start\", or \"Restart TDDFT\". Terminating.\n");

    If.RegisterInputKey("freeze_orbital_step", &lc.freeze_orbital_step, 0, 100000, 90, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "freeze orbital after this step", 
                     "");
    If.RegisterInputKey("freeze_rho_steps", &lc.freeze_rho_steps, 0, 100, 0, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "freeze rho which read from RMG for a number of steps",""); 

    If.RegisterInputKey("orbital_mixing_method", NULL, &lc.orbital_mixing_method, "Pulay",
                     CHECK_AND_TERMINATE, REQUIRED, mg_method,
                     "mixing type for orbitals  \"Steepest Descent\", \"Pulay\", or \"KAIN\".\n", 
                     "start_mode must be one of \"Steepest Descent\", \"Pulay\", or \"KAIN\". Terminating.\n");

    If.RegisterInputKey("orbital_pulay_order", &lc.orbital_pulay_order, 0, 10, 2, 
                     CHECK_AND_FIX, OPTIONAL, "", "");

    If.RegisterInputKey("orbital_pulay_refresh", &lc.orbital_pulay_refresh, 0, 100, 100, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("orbital_mixing", &lc.orbital_pulay_mixfirst, 0.0, 1.0, 0.1,
                     CHECK_AND_FIX, OPTIONAL,
                     "mixing parameter when linear mixing is used or first step in Pulay mixing.\n",
                     "must lie in the range (0.0, 1.0) Resetting to the default value of 0.1.\n");
    If.RegisterInputKey("orbital_pulay_scale", &lc.orbital_pulay_scale, 0.0, 1.0, 0.1,
                     CHECK_AND_FIX, OPTIONAL,
                     "mix orbitals and their residuals in Pulay mixing.\n",
                     "must lie in the range (0.0, 1.0) Resetting to the default value of 0.1.\n");


    If.RegisterInputKey("do_movable_orbital_centers", &lc.movingCenter, false, "");
    If.RegisterInputKey("band_width_reduction", &lc.bandwidthreduction, false, "");

    If.RegisterInputKey("movable_orbital_centers_steps", &lc.movingSteps, 1, 10000, 40, 
                     CHECK_AND_FIX, OPTIONAL, "", "");

    If.RegisterInputKey("number_of_orbitals", &lc.num_states, 1, INT_MAX, 4, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("number_of_atoms", &lc.num_ions, 1, INT_MAX, 1, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("atomic_orbital_files", NULL, "", CHECK_AND_FIX, OPTIONAL, "","");
    If.RegisterInputKey("ReadFromRMG", &lc.ON_read_from_RMG, false, "");

    If.RegisterInputKey("InputFileFromRMG", &RMGfile, "RMG/Waves/wave.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Output file/path to store wavefunctions and other binary data.\n",
                     "");



    If.LoadInputKeys();

    if( (!RMGfile.length() ) && lc.ON_read_from_RMG ) 
    {
        printf(" \nNeed a file name and path from RMG run, so we can read potentials and charge density as initial ones\n");
        fflush(NULL);
        exit(0);
    }
    else
    {
        std::strncpy(lc.infile_ON_from_RMG, RMGfile.c_str(), sizeof(lc.infile_ON_from_RMG));
    }

    InputKey *Key;
    std::vector<std::string> lines;
    std::string delims = "\r\n^";
    std::string field_delims = " \t";
    Key = NewMap.at("atomic_orbital_files");
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
                if(!std::strcmp(fields[0].c_str(), Species[isp].atomic_symbol)) {
                    lc.file_atomic_orbit.push_back(fields[1]);
                }
            }

        }

    }

}

