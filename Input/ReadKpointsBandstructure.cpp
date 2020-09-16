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
    std::vector<std::string> KpointList;
    std::unordered_map<std::string, InputKey *> NewMap;
    int nkpts;

    RmgInputFile If(cfile, NewMap, pct.img_comm);

    If.RegisterInputKey("kpoints_bandstructure", &KpointArray, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "List of kpoints to use in a bandstructure calculation. For more detailed information look at the github wiki page on kpoint calculations.\n", "");
    
    
    If.LoadInputKeys();
    // Process atoms
    boost::trim(KpointArray);
    boost::trim_if(KpointArray, boost::algorithm::is_any_of("\"^"));

    boost::algorithm::split( KpointList, KpointArray, boost::is_any_of(line_delims), boost::token_compress_on );

    
    std::vector<std::string>::iterator it, it1;
    int num_lines = KpointList.size();
    if(num_lines < 2) return 0;
    nkpts=0;

    std::string Kpoint ;
    std::vector<std::string> KpointComponents;
    for (int kpt = 1; kpt < num_lines; kpt++)
    {

        Kpoint = KpointList[kpt];
        boost::trim(Kpoint);

        boost::algorithm::split( KpointComponents, Kpoint, boost::is_any_of(whitespace_delims), boost::token_compress_on );

        size_t ncomp = KpointComponents.size();
        if((ncomp != 5) ) throw RmgFatalException() << "Synax error in kpoint_band information near " << Kpoint << "\n";

        nkpts += std::atoi(KpointComponents[3].c_str());

    }

    lc.num_kpts = nkpts + 1;
    lc.kp.resize(lc.num_kpts);

    for(int kpt=0; kpt < lc.num_kpts; kpt++) 
        strcpy(ct.kp[kpt].symbol,  "");

    Kpoint = KpointList[0];
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

        std::string Kpoint = KpointList[kpt];
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

    for (int kpt = 0; kpt < ct.num_kpts; kpt++) {
        double v1, v2, v3;
        v1 = 0.0;
        v2 = 0.0;
        v3 = 0.0;

        for(int ir = 0; ir<3; ir++)
        {
            v1 = ct.kp[kpt].kpt[0] *Rmg_L.b0[0]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[0] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[0];
            v2 = ct.kp[kpt].kpt[0] *Rmg_L.b0[1]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[1] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[1];
            v3 = ct.kp[kpt].kpt[0] *Rmg_L.b0[2]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[2] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[2];
        }

        ct.kp[kpt].kvec[0] = v1 * twoPI;
        ct.kp[kpt].kvec[1] = v2 * twoPI;
        ct.kp[kpt].kvec[2] = v3 * twoPI;
        ct.kp[kpt].kmag = (v1 * v1 + v2 * v2 + v3 * v3) * twoPI * twoPI;

    }


    return 1;

//    rmg_printf("\n num_k %d", ct.num_kpts);
//    for(int kpt = 0; kpt < ct.num_kpts; kpt++)
//       rmg_printf("\n kvec %d  %f %f %f %s", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].symbol);
//    rmg_printf("\n");

//    exit(0);

}
