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



/**********************************************************************

    Column 1: kx-coordinate
    Column 2: ky-coordinate
    Column 3: kz-coordinate
    Column 4: num of kpoint inserted between current point and previous point
    Column 5: kpoint symbol for plot
            

Example:

kpoints = "
    0.0  0.0  0.0   0   \xG
    0.5  0.5  0.5   21  L
"


    cfile        = name of the file containing the kpoint info
    InputMap     = Control Map. May not be needed by all atomic input
                   drivers but is useful for reading the RMG format.
    
**********************************************************************/

namespace Ri = RmgInput;

int ReadKpointsBandstructure(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::string KpointArray;
    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::vector<std::string> Kpoints;
    std::unordered_map<std::string, InputKey *> NewMap;
    int nkpts;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("kpoints_bandstructure", &KpointArray, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "kpoints list \n",
                     "");
    
    
    If.LoadInputKeys();
    // Process atoms
    boost::trim(KpointArray);
    boost::trim_if(KpointArray, boost::algorithm::is_any_of("\"^"));

    boost::algorithm::split( Kpoints, KpointArray, boost::is_any_of(line_delims), boost::token_compress_on );

    
    std::vector<std::string>::iterator it, it1;
    int num_lines = Kpoints.size();
    if(num_lines < 2) return 0;
    nkpts=0;

    std::string Kpoint ;
    std::vector<std::string> KpointComponents;
    for (int kpt = 1; kpt < num_lines; kpt++)
    {

        Kpoint = Kpoints[kpt];
        boost::trim(Kpoint);

        boost::algorithm::split( KpointComponents, Kpoint, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = KpointComponents.size();
        if((ncomp != 5) ) throw RmgFatalException() << "Synax error in kpoint_band information near " << Kpoint << "\n";

        nkpts += std::atoi(KpointComponents[3].c_str());

    }

    lc.num_kpts = nkpts + 1;
    lc.kp = new KPOINT[lc.num_kpts]();

    for(int kpt=0; kpt < lc.num_kpts; kpt++) 
        strcpy(ct.kp[kpt].symbol,  "");

    Kpoint = Kpoints[0];
    boost::trim(Kpoint);

    boost::algorithm::split( KpointComponents, Kpoint, boost::is_any_of(whitespace_delims), boost::token_compress_on );

    double kx0, ky0, kz0, kx1, ky1, kz1, dx, dy, dz;
    kx0 = std::atof(KpointComponents[0].c_str());
    ky0 = std::atof(KpointComponents[1].c_str());
    kz0 = std::atof(KpointComponents[2].c_str());


    nkpts=0;
    strcpy(ct.kp[nkpts].symbol,  KpointComponents[4].c_str());
    for (int kpt = 1; kpt < num_lines; kpt++)
    {

        std::string Kpoint = Kpoints[kpt];
        boost::trim(Kpoint);

        std::vector<std::string> KpointComponents;
        boost::algorithm::split( KpointComponents, Kpoint, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        kx1 = std::atof(KpointComponents[0].c_str());
        ky1 = std::atof(KpointComponents[1].c_str());
        kz1 = std::atof(KpointComponents[2].c_str());
        int num = std::atoi(KpointComponents[3].c_str());
        dx = (kx1 - kx0)/ num; 
        dy = (ky1 - ky0)/ num; 
        dz = (kz1 - kz0)/ num; 
        for(int i = 0; i < num; i++)
        {
            ct.kp[nkpts].kpt[0] = kx0 + dx * i;
            ct.kp[nkpts].kpt[1] = ky0 + dy * i;
            ct.kp[nkpts].kpt[2] = kz0 + dz * i;
            nkpts++;
        }
        strcpy(ct.kp[nkpts].symbol,  KpointComponents[4].c_str());
        kx0 = kx1;
        ky0 = ky1;
        kz0 = kz1;

    }
    ct.kp[nkpts].kpt[0] = kx0;
    ct.kp[nkpts].kpt[1] = ky0;
    ct.kp[nkpts].kpt[2] = kz0;
    return 1;

//    rmg_printf("\n num_k %d", ct.num_kpts);
//    for(int kpt = 0; kpt < ct.num_kpts; kpt++)
//       rmg_printf("\n kvec %d  %f %f %f %s", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].symbol);
//    rmg_printf("\n");


//    exit(0);

}



