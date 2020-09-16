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

    Reads internal RMG atomic format which consists of a single atom
    per line with the entire block of lines enclosed in double quotes
    and the block identified by an atoms keyword.

    Column 1: kx-coordinate
    Column 2: ky-coordinate
    Column 3: kz-coordinate
    Column 4: weight
            

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

void ReadKpoints(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    std::string KpointArray;
    std::string line_delims = "^\n";
    std::string whitespace_delims = " \n\t";
    std::vector<std::string> KpointList;
    std::unordered_map<std::string, InputKey *> NewMap;
    int nkpts;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("kpoints", &KpointArray, "",
                     CHECK_AND_FIX, REQUIRED,
                     "Normally kpoints are specified using the kpoint_mesh and kpoint_is_shift options but one can also enter a list of kpoints and their weights with this option. If kpoint_mesh is not specified or this is a bandstructure calculation this is required otherwise it is optional. \n",
                     "");
    
    
    If.LoadInputKeys();
    // Process atoms
    boost::trim(KpointArray);
    boost::trim_if(KpointArray, boost::algorithm::is_any_of("\"^"));

    boost::algorithm::split( KpointList, KpointArray, boost::is_any_of(line_delims), boost::token_compress_on );

    


    lc.num_kpts = KpointList.size();
    lc.kp.resize(lc.num_kpts);

    std::vector<std::string>::iterator it, it1;
    nkpts=0;

    for (it = KpointList.begin(); it != KpointList.end(); ++it) {

        std::string Kpoint = *it;
        boost::trim(Kpoint);

        std::vector<std::string> KpointComponents;
        boost::algorithm::split( KpointComponents, Kpoint, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = KpointComponents.size();
        if((ncomp != 4) ) throw RmgFatalException() << "Synax error in kpoint information near " << Kpoint << "\n";

        // First field should be an atomic symbol
        it1 = KpointComponents.begin();
      
        std::string xstr = *it1;
        lc.kp[nkpts].kpt[0] = std::atof(xstr.c_str());
        it1++;
        std::string ystr = *it1;
        lc.kp[nkpts].kpt[1] = std::atof(ystr.c_str());
        it1++;
        std::string zstr = *it1;
        lc.kp[nkpts].kpt[2] = std::atof(zstr.c_str());
        it1++;
        std::string wstr = *it1;
        lc.kp[nkpts].kweight = std::atof(wstr.c_str());

        nkpts++;

    }

    lc.is_gamma = true;
    for (int kpt = 0; kpt < lc.num_kpts; kpt++) {
        double v1, v2, v3;
        v1 = 0.0;
        v2 = 0.0;
        v3 = 0.0;

        for(int ir = 0; ir<3; ir++)
        {
            v1 = lc.kp[kpt].kpt[0] *Rmg_L.b0[0]
                + lc.kp[kpt].kpt[1] *Rmg_L.b1[0] 
                + lc.kp[kpt].kpt[2] *Rmg_L.b2[0];
            v2 = lc.kp[kpt].kpt[0] *Rmg_L.b0[1]
                + lc.kp[kpt].kpt[1] *Rmg_L.b1[1] 
                + lc.kp[kpt].kpt[2] *Rmg_L.b2[1];
            v3 = lc.kp[kpt].kpt[0] *Rmg_L.b0[2]
                + lc.kp[kpt].kpt[1] *Rmg_L.b1[2] 
                + lc.kp[kpt].kpt[2] *Rmg_L.b2[2];
        }

        lc.kp[kpt].kvec[0] = v1 * twoPI;
        lc.kp[kpt].kvec[1] = v2 * twoPI;
        lc.kp[kpt].kvec[2] = v3 * twoPI;
        lc.kp[kpt].kmag = (v1 * v1 + v2 * v2 + v3 * v3) * twoPI * twoPI;

        if(lc.kp[kpt].kmag != 0.0) lc.is_gamma = false;
    }
        
    //rmg_printf("\n num_k %d", ct.num_kpts);
    //for(int kpt = 0; kpt < ct.num_kpts; kpt++)
     //   rmg_printf("\n kvec %d  %f %f %f %f", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);
    //rmg_printf("\n");


}

