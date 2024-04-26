#if INTERNAL_PP

#include <unordered_map>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp> 
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/write.hpp>
#include <boost/iostreams/concepts.hpp>   
#include <boost/iostreams/device/back_inserter.hpp>
#include <iterator> 
#include <initializer_list>
#include "InternalPseudo.h"
#include "pseudo_list_ncpp_standard.h"
#include "rmg_error.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "main.h"

std::unordered_map<std::string, compressed_pp> NCPP_STANDARD_FILES ({
{"AG"    , Ag_standard_upf_bz2},
{"AL"    , Al_standard_upf_bz2},
{"AR"    , Ar_standard_upf_bz2},
{"AS"    , As_standard_upf_bz2},
{"AU"    , Au_standard_upf_bz2},
{"BA"    , Ba_standard_upf_bz2},
{"BE"    , Be_standard_upf_bz2},
{"BI"    , Bi_standard_upf_bz2},
{"BR"    , Br_standard_upf_bz2},
{"B"    , B_standard_upf_bz2},
{"CA"    , Ca_standard_upf_bz2},
{"CD"    , Cd_standard_upf_bz2},
{"CL"    , Cl_standard_upf_bz2},
{"CO"    , Co_standard_upf_bz2},
{"CR"    , Cr_standard_upf_bz2},
{"CS"    , Cs_standard_upf_bz2},
{"C"    , C_standard_upf_bz2},
{"CU"    , Cu_standard_upf_bz2},
{"FE"    , Fe_standard_upf_bz2},
{"F"    , F_standard_upf_bz2},
{"GA"    , Ga_standard_upf_bz2},
{"GE"    , Ge_standard_upf_bz2},
{"HE"    , He_standard_upf_bz2},
{"HF"    , Hf_standard_upf_bz2},
{"HG"    , Hg_standard_upf_bz2},
{"H"    , H_standard_upf_bz2},
{"IN"    , In_standard_upf_bz2},
{"IR"    , Ir_standard_upf_bz2},
{"I"    , I_standard_upf_bz2},
{"KR"    , Kr_standard_upf_bz2},
{"K"    , K_standard_upf_bz2},
{"LA"    , La_standard_upf_bz2},
{"LI"    , Li_standard_upf_bz2},
{"LU"    , Lu_standard_upf_bz2},
{"MG"    , Mg_standard_upf_bz2},
{"MN"    , Mn_standard_upf_bz2},
{"MO"    , Mo_standard_upf_bz2},
{"NA"    , Na_standard_upf_bz2},
{"NB"    , Nb_standard_upf_bz2},
{"NE"    , Ne_standard_upf_bz2},
{"NI"    , Ni_standard_upf_bz2},
{"N"    , N_standard_upf_bz2},
{"OS"    , Os_standard_upf_bz2},
{"O"    , O_standard_upf_bz2},
{"PB"    , Pb_standard_upf_bz2},
{"PD"    , Pd_standard_upf_bz2},
{"PO"    , Po_standard_upf_bz2},
{"P"    , P_standard_upf_bz2},
{"PT"    , Pt_standard_upf_bz2},
{"RB"    , Rb_standard_upf_bz2},
{"RE"    , Re_standard_upf_bz2},
{"RH"    , Rh_standard_upf_bz2},
{"RN"    , Rn_standard_upf_bz2},
{"RU"    , Ru_standard_upf_bz2},
{"SB"    , Sb_standard_upf_bz2},
{"SC"    , Sc_standard_upf_bz2},
{"SE"    , Se_standard_upf_bz2},
{"SI"    , Si_standard_upf_bz2},
{"SN"    , Sn_standard_upf_bz2},
{"SR"    , Sr_standard_upf_bz2},
{"S"    , S_standard_upf_bz2},
{"TA"    , Ta_standard_upf_bz2},
{"TC"    , Tc_standard_upf_bz2},
{"TE"    , Te_standard_upf_bz2},
{"TI"    , Ti_standard_upf_bz2},
{"TL"    , Tl_standard_upf_bz2},
{"V"    , V_standard_upf_bz2},
{"W"    , W_standard_upf_bz2},
{"XE"    , Xe_standard_upf_bz2},
{"Y"    , Y_standard_upf_bz2},
{"ZN"    , Zn_standard_upf_bz2},
{"ZR"    , Zr_standard_upf_bz2}
});

std::unordered_map<std::string, unsigned int> NCPP_STANDARD_FILES_LEN ({
{"AG" , 69238},
{"AL" , 60998},
{"AR" , 34074},
{"AS" , 54324},
{"AU" , 64013},
{"BA" , 78114},
{"BE" , 41242},
{"BI" , 62654},
{"BR" , 39874},
{"B" , 47101},
{"CA" , 65413},
{"CD" , 61732},
{"CL" , 36701},
{"CO" , 62713},
{"CR" , 69282},
{"CS" , 94410},
{"C" , 39245},
{"CU" , 66343},
{"FE" , 62640},
{"F" , 29475},
{"GA" , 72074},
{"GE" , 61389},
{"HE" , 18691},
{"HF" , 78991},
{"HG" , 60294},
{"H" , 26814},
{"IN" , 76082},
{"IR" , 64244},
{"I" , 47323},
{"KR" , 37210},
{"K" , 75761},
{"LA" , 84442},
{"LI" , 53282},
{"LU" , 100243},
{"MG" , 53212},
{"MN" , 64633},
{"MO" , 71954},
{"NA" , 63228},
{"NB" , 74533},
{"NE" , 25570},
{"NI" , 62052},
{"N" , 34048},
{"OS" , 65546},
{"O" , 31878},
{"PB" , 69197},
{"PD" , 49949},
{"PO" , 57163},
{"P" , 44130},
{"PT" , 65944},
{"RB" , 82863},
{"RE" , 67041},
{"RH" , 69111},
{"RN" , 51709},
{"RU" , 68590},
{"SB" , 57868},
{"SC" , 70780},
{"SE" , 49580},
{"SI" , 51022},
{"SN" , 65257},
{"SR" , 69320},
{"S" , 39416},
{"TA" , 70348},
{"TC" , 66749},
{"TE" , 56002},
{"TI" , 68709},
{"TL" , 79998},
{"V" , 66760},
{"W" , 67659},
{"XE" , 43721},
{"Y" , 78054},
{"ZN" , 61288},
{"ZR" , 73051}
});
namespace io = boost::iostreams;

std::string GetInternalPseudo_ncpp_standard(const char *symbol)
{

    std::string decompressed;
    std::string lsym(symbol);
    boost::to_upper(lsym);

    char *pptr = (char *)NCPP_STANDARD_FILES[lsym];
    unsigned int len = NCPP_STANDARD_FILES_LEN[lsym];
    std::vector<char> pp_file_v(pptr, pptr + len);

    io::filtering_ostream os;
    os.push(io::bzip2_decompressor());
    os.push(io::back_inserter(decompressed));

    os.flush();
    io::write(os, pptr, (std::streamsize)pp_file_v.size());
    os.flush();
    io::close(os);
 
    std::string pp_file = lsym + std::string("_standard_rmg_internal.upf");
    if(pct.worldrank == 0)
    {
        int fhand;
        // If we can't write this what should we do?
        fhand = open(pp_file.c_str(), O_RDWR|O_CREAT|O_TRUNC, S_IREAD | S_IWRITE);
        if(fhand < 0)
            rmg_error_handler (__FILE__, __LINE__, " Error saving pseudopotential file. Terminating.");

        write(fhand, decompressed.c_str(), decompressed.length());
        close(fhand);
    }
    return decompressed;

    
}

#endif
