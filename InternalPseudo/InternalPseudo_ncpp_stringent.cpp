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
#include "pseudo_list_ncpp_stringent.h"

std::unordered_map<std::string, compressed_pp> NCPP_STRINGENT_FILES ({
{"AG"   , Ag_stringent_upf_bz2},
{"AL"   , Al_stringent_upf_bz2},
{"AR"   , Ar_stringent_upf_bz2},
{"AS"   , As_stringent_upf_bz2},
{"AU"   , Au_stringent_upf_bz2},
{"BA"   , Ba_stringent_upf_bz2},
{"BE"   , Be_stringent_upf_bz2},
{"BI"   , Bi_stringent_upf_bz2},
{"BR"   , Br_stringent_upf_bz2},
{"B"   , B_stringent_upf_bz2},
{"CA"   , Ca_stringent_upf_bz2},
{"CD"   , Cd_stringent_upf_bz2},
{"CL"   , Cl_stringent_upf_bz2},
{"CO"   , Co_stringent_upf_bz2},
{"CR"   , Cr_stringent_upf_bz2},
{"CS"   , Cs_stringent_upf_bz2},
{"C"   , C_stringent_upf_bz2},
{"CU"   , Cu_stringent_upf_bz2},
{"FE"   , Fe_stringent_upf_bz2},
{"F"   , F_stringent_upf_bz2},
{"GA"   , Ga_stringent_upf_bz2},
{"GE"   , Ge_stringent_upf_bz2},
{"HE"   , He_stringent_upf_bz2},
{"HF"   , Hf_stringent_upf_bz2},
{"HG"   , Hg_stringent_upf_bz2},
{"H"   , H_stringent_upf_bz2},
{"IN"   , In_stringent_upf_bz2},
{"IR"   , Ir_stringent_upf_bz2},
{"I"   , I_stringent_upf_bz2},
{"KR"   , Kr_stringent_upf_bz2},
{"K"   , K_stringent_upf_bz2},
{"LA"   , La_stringent_upf_bz2},
{"LI"   , Li_stringent_upf_bz2},
{"LU"   , Lu_stringent_upf_bz2},
{"MG"   , Mg_stringent_upf_bz2},
{"MN"   , Mn_stringent_upf_bz2},
{"MO"   , Mo_stringent_upf_bz2},
{"NA"   , Na_stringent_upf_bz2},
{"NB"   , Nb_stringent_upf_bz2},
{"NE"   , Ne_stringent_upf_bz2},
{"NI"   , Ni_stringent_upf_bz2},
{"N"   , N_stringent_upf_bz2},
{"OS"   , Os_stringent_upf_bz2},
{"O"   , O_stringent_upf_bz2},
{"PB"   , Pb_stringent_upf_bz2},
{"PD"   , Pd_stringent_upf_bz2},
{"PO"   , Po_stringent_upf_bz2},
{"P"   , P_stringent_upf_bz2},
{"PT"   , Pt_stringent_upf_bz2},
{"RB"   , Rb_stringent_upf_bz2},
{"RE"   , Re_stringent_upf_bz2},
{"RH"   , Rh_stringent_upf_bz2},
{"RN"   , Rn_stringent_upf_bz2},
{"RU"   , Ru_stringent_upf_bz2},
{"SB"   , Sb_stringent_upf_bz2},
{"SC"   , Sc_stringent_upf_bz2},
{"SE"   , Se_stringent_upf_bz2},
{"SI"   , Si_stringent_upf_bz2},
{"SN"   , Sn_stringent_upf_bz2},
{"SR"   , Sr_stringent_upf_bz2},
{"S"   , S_stringent_upf_bz2},
{"TA"   , Ta_stringent_upf_bz2},
{"TC"   , Tc_stringent_upf_bz2},
{"TE"   , Te_stringent_upf_bz2},
{"TI"   , Ti_stringent_upf_bz2},
{"TL"   , Tl_stringent_upf_bz2},
{"V"   , V_stringent_upf_bz2},
{"U"   , U_ONCV_PBE_sr_upf_bz2},
{"W"   , W_stringent_upf_bz2},
{"XE"   , Xe_stringent_upf_bz2},
{"Y"   , Y_stringent_upf_bz2},
{"ZN"   , Zn_stringent_upf_bz2},
{"ZR"   , Zr_stringent_upf_bz2}
});

std::unordered_map<std::string, unsigned int> NCPP_STRINGENT_FILES_LEN ({
{"AG" ,  69238},
{"AL" ,  60920},
{"AR" ,  34074},
{"AS" ,  69768},
{"AU" ,  64013},
{"BA" ,  75745},
{"BE" ,  40372},
{"BI" ,  77165},
{"BR" ,  46048},
{"B" ,  47101},
{"CA" ,  65413},
{"CD" ,  60492},
{"CL" ,  36699},
{"CO" ,  61145},
{"CR" ,  65142},
{"CS" ,  94410},
{"C" ,  39245},
{"CU" ,  64283},
{"FE" ,  60318},
{"F" ,  29475},
{"GA" ,  87691},
{"GE" ,  76358},
{"HE" ,  18691},
{"HF" ,  86602},
{"HG" ,  59885},
{"H" ,  26814},
{"IN" ,  94795},
{"IR" ,  64244},
{"I" ,  49609},
{"KR" ,  37210},
{"K" ,  75761},
{"LA" ,  84442},
{"LI" ,  51991},
{"LU" ,  100243},
{"MG" ,  53212},
{"MN" ,  61040},
{"MO" ,  71954},
{"NA" ,  63228},
{"NB" ,  74533},
{"NE" ,  25305},
{"NI" ,  61134},
{"N" ,  34048},
{"OS" ,  65546},
{"O" ,  31622},
{"PB" ,  83814},
{"PD" ,  49949},
{"PO" ,  70166},
{"P" ,  43954},
{"PT" ,  65944},
{"RB" ,  83105},
{"RE" ,  67041},
{"RH" ,  69111},
{"RN" ,  48922},
{"RU" ,  68590},
{"SB" ,  73471},
{"SC" ,  70780},
{"SE" ,  63834},
{"SI" ,  51027},
{"SN" ,  79616},
{"SR" ,  69163},
{"S" ,  39817},
{"TA" ,  78972},
{"TC" ,  66749},
{"TE" ,  67834},
{"TI" ,  68709},
{"TL" ,  93279},
{"U" , 414920},
{"V" ,  66760},
{"W" ,  67659},
{"XE" ,  41710},
{"Y" ,  78054},
{"ZN" ,  61398},
{"ZR" ,  73051}
});
namespace io = boost::iostreams;

std::string GetInternalPseudo_ncpp_stringent(const char *symbol)
{

    std::string decompressed;
    std::string lsym(symbol);
    boost::to_upper(lsym);

    char *pptr = (char *)NCPP_STRINGENT_FILES[lsym];
    unsigned int len = NCPP_STRINGENT_FILES_LEN[lsym];
    std::vector<char> pp_file_v(pptr, pptr + len);

    io::filtering_ostream os;
    os.push(io::bzip2_decompressor());
    os.push(io::back_inserter(decompressed));

    io::write(os, pptr, (std::streamsize)pp_file_v.size());
    return decompressed;

    
}

#endif
