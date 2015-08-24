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
#include "LCR.h"
#include "twoParts.h"



/**********************************************************************

    read other stuff for ON branch

    cfile        = name of the file containing the kpoint info
    InputMap     = Control Map. May not be needed by all atomic input
                   drivers but is useful for reading the RMG format.
    
**********************************************************************/

namespace Ri = RmgInput;

void ReadBranchNEGF(char *cfile, CONTROL& lc, complex_energy_integral& cei, COMPASS& potcompass, COMPASS& rhocompass)
{

    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::unordered_map<std::string, InputKey *> NewMap;

    int negf_runflag;

    std::string BlockDim; 
    std::string PotCompass;
    std::string RhoCompass;
    std::string AveragePlane;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    


    If.RegisterInputKey("start_mode", NULL, &lc.runflag, "",
                     CHECK_AND_TERMINATE, REQUIRED, start_mode,
                     "Type of run. Choices are \"Random Start\", \"Restart From File\", \"FIREBALL Start\", \"Gaussian Start\", or \"Restart TDDFT\".\n", 
                     "start_mode must be one of  \"Random Start\", \"Restart From File\", \"FIREBALL Start\", \"Gaussian Start\", or \"Restart TDDFT\". Terminating.\n");


    If.RegisterInputKey("number_of_orbitals", &lc.num_states, 1, INT_MAX, 4, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("number_of_atoms", &lc.num_ions, 1, INT_MAX, 1, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("energy_point_insert", NULL, &cei.energy_point_insert, "None",
                     CHECK_AND_TERMINATE, OPTIONAL, energy_point_insert_mode,
                     "Type of Choices are \"None\", \"Simpson\", or \"Sharp Peaks\".\n", 
                     "energy point insert must be one of \"None\", \"Simpson\", or \"Sharp Peaks\". Terminating.\n");
    If.RegisterInputKey("Simpson_depth", &lc.simpson_depth, 0, 10, 0, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("Simpson_tol", &lc.simpson_tol, 1e-10, 1.0, 1e-3, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("probe_noneq", &cei.probe_noneq, 1, 100, 1, 
                     CHECK_AND_FIX, OPTIONAL, "", "");

    If.RegisterInputKey("gbias_begin", &lc.gbias_begin, -DBL_MAX, DBL_MAX, -100.0,
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("gbias_end", &lc.gbias_end, -DBL_MAX, DBL_MAX, 622.0,
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("BT", &ct.BT, -DBL_MAX, DBL_MAX, 4.0,
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("gate_bias", &lc.gate_bias, -DBL_MAX, DBL_MAX, 0.0,
                     CHECK_AND_FIX, OPTIONAL, "", "");

    If.RegisterInputKey("vcomp_Lbegin", &lc.vcomp_Lbegin, -INT_MAX, INT_MAX, -1, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("vcomp_Lend", &lc.vcomp_Lend, -INT_MAX, INT_MAX, -1, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("vcomp_Rbegin", &lc.vcomp_Rbegin, -INT_MAX, INT_MAX, -1, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("vcomp_Rend", &lc.vcomp_Rend, -INT_MAX, INT_MAX, -1, 
                     CHECK_AND_FIX, OPTIONAL, "", "");

    If.RegisterInputKey("auto_3Ddos", &lc.auto_3Ddos, -INT_MAX, INT_MAX, 0, 
                     CHECK_AND_FIX, OPTIONAL, "", "");

    If.RegisterInputKey("metalic", &lc.metal, true, "");

    If.RegisterInputKey("num_blocks", &lc.num_blocks, 3, INT_MAX, 3, 
                     CHECK_AND_FIX, OPTIONAL, "", "");
    If.RegisterInputKey("blocks_dim", &BlockDim, "", CHECK_AND_FIX, REQUIRED,"","");
    If.RegisterInputKey("potential_compass", &PotCompass, "", CHECK_AND_FIX, REQUIRED, "","");
    If.RegisterInputKey("chargedensity_compass", &RhoCompass, "",  CHECK_AND_FIX, REQUIRED, "","");

    If.RegisterInputKey("start_mode_NEGF", &negf_runflag, 0, INT_MAX, 112, 
                     CHECK_AND_FIX, REQUIRED, "", "");
    If.RegisterInputKey("average_plane_rho", &AveragePlane, "1 0 0 0 1", CHECK_AND_FIX, OPTIONAL, "","");


    If.LoadInputKeys();

    lc.runflag = negf_runflag;
    char *tbuf = new char[255];
    strcpy(tbuf, PotCompass.c_str());
    potcompass.type = std::strtol(tbuf, &tbuf, 10);
    potcompass.box1.x1 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    potcompass.box1.x2 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    potcompass.box1.y1 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    potcompass.box1.y2 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    potcompass.box1.z1 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    potcompass.box1.z2 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();

    strcpy(tbuf, RhoCompass.c_str());
    rhocompass.type = std::strtol(tbuf, &tbuf, 10);
    rhocompass.box1.x1 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    rhocompass.box1.x2 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    rhocompass.box1.y1 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    rhocompass.box1.y2 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    rhocompass.box1.z1 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();
    rhocompass.box1.z2 = std::strtol(tbuf, &tbuf, 10) * get_FG_RATIO();

    strcpy(tbuf, AveragePlane.c_str());
    for(int i = 0; i < 5; i++)
        ct.plane[i] = std::strtol(tbuf, &tbuf, 10);

    strcpy(tbuf, BlockDim.c_str());
    int n, num_st = 0, n_block = 0;
    while( (n = std::strtol(tbuf, &tbuf, 10)) > 0)
    {  
        ct.block_dim[n_block] = n;
        n_block++;
        num_st += n;
    }

    if((n_block != lc.num_blocks) | (num_st != ct.num_states) )
       throw RmgFatalException()<<"block_dim is not right\n"<< BlockDim.c_str()<<"\n";


}
