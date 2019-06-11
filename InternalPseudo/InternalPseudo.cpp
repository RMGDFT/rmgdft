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
#include "pseudo_list.h"



std::unordered_map<std::string, compressed_pp> USPP_FILES ({
{"AG"   , ag_pbe_v1_4_uspp_F_UPF_bz2},
{"AL"   , al_pbe_v1_uspp_F_UPF_bz2},
{"AS"   , as_pbe_v1_uspp_F_UPF_bz2},
{"AU"   , au_pbe_v1_uspp_F_UPF_bz2},
{"BA"   , ba_pbe_v1_uspp_F_UPF_bz2},
{"BE"   , be_pbe_v1_4_uspp_F_UPF_bz2},
{"B"    , b_pbe_v1_4_uspp_F_UPF_bz2},
{"BR"   , br_pbe_v1_4_uspp_F_UPF_bz2},
{"CA"   , ca_pbe_v1_uspp_F_UPF_bz2},
{"CD"   , cd_pbe_v1_uspp_F_UPF_bz2},
{"CL"   , cl_pbe_v1_4_uspp_F_UPF_bz2},
{"CO"   , co_pbe_v1_2_uspp_F_UPF_bz2},
{"C"    , c_pbe_v1_2_uspp_F_UPF_bz2},
{"CR"   , cr_pbe_v1_5_uspp_F_UPF_bz2},
{"CS"   , cs_pbe_v1_uspp_F_UPF_bz2},
{"CU"   , cu_pbe_v1_2_uspp_F_UPF_bz2},
{"FE"   , fe_pbe_v1_5_uspp_F_UPF_bz2},
{"F"    , f_pbe_v1_4_uspp_F_UPF_bz2},
{"GA"   , ga_pbe_v1_4_uspp_F_UPF_bz2},
{"GE"   , ge_pbe_v1_4_uspp_F_UPF_bz2},
{"HF"   , hf_pbe_v1_uspp_F_UPF_bz2},
{"HG"   , hg_pbe_v1_uspp_F_UPF_bz2},
{"H"    , h_pbe_v1_4_uspp_F_UPF_bz2},
{"IN"   , in_pbe_v1_4_uspp_F_UPF_bz2},
{"I"    , i_pbe_v1_uspp_F_UPF_bz2},
{"IR"   , ir_pbe_v1_2_uspp_F_UPF_bz2},
{"K"    , k_pbe_v1_4_uspp_F_UPF_bz2},
{"LA"   , la_pbe_v1_uspp_F_UPF_bz2},
{"LI"   , li_pbe_v1_4_uspp_F_UPF_bz2},
{"MG"   , mg_pbe_v1_4_uspp_F_UPF_bz2},
{"MN"   , mn_pbe_v1_5_uspp_F_UPF_bz2},
{"MO"   , mo_pbe_v1_uspp_F_UPF_bz2},
{"NA"   , na_pbe_v1_5_uspp_F_UPF_bz2},
{"NB"   , nb_pbe_v1_uspp_F_UPF_bz2},
{"NI"   , ni_pbe_v1_4_uspp_F_UPF_bz2},
{"N"    , n_pbe_v1_2_uspp_F_UPF_bz2},
{"O"    , o_pbe_v1_2_uspp_F_UPF_bz2},
{"OS"   , os_pbe_v1_2_uspp_F_UPF_bz2},
{"PB"   , pb_pbe_v1_uspp_F_UPF_bz2},
{"PD"   , pd_pbe_v1_4_uspp_F_UPF_bz2},
{"P"    , p_pbe_v1_5_uspp_F_UPF_bz2},
{"PT"   , pt_pbe_v1_4_uspp_F_UPF_bz2},
{"RB"   , rb_pbe_v1_uspp_F_UPF_bz2},
{"RE"   , re_pbe_v1_2_uspp_F_UPF_bz2},
{"RH"   , rh_pbe_v1_4_uspp_F_UPF_bz2},
{"RU"   , ru_pbe_v1_2_uspp_F_UPF_bz2},
{"SB"   , sb_pbe_v1_4_uspp_F_UPF_bz2},
{"SC"   , sc_pbe_v1_uspp_F_UPF_bz2},
{"SE"   , se_pbe_v1_uspp_F_UPF_bz2},
{"SI"   , si_pbe_v1_uspp_F_UPF_bz2},
{"SN"   , sn_pbe_v1_4_uspp_F_UPF_bz2},
{"S"    , s_pbe_v1_4_uspp_F_UPF_bz2},
{"SR"   , sr_pbe_v1_uspp_F_UPF_bz2},
{"TA"   , ta_pbe_v1_uspp_F_UPF_bz2},
{"TC"   , tc_pbe_v1_uspp_F_UPF_bz2},
{"TE"   , te_pbe_v1_uspp_F_UPF_bz2},
{"TI"   , ti_pbe_v1_4_uspp_F_UPF_bz2},
{"TL"   , tl_pbe_v1_2_uspp_F_UPF_bz2},
{"V"    , v_pbe_v1_4_uspp_F_UPF_bz2},
{"W"    , w_pbe_v1_2_uspp_F_UPF_bz2},
{"Y"    , y_pbe_v1_4_uspp_F_UPF_bz2},
{"ZN"   , zn_pbe_v1_uspp_F_UPF_bz2},
{"ZR"   , zr_pbe_v1_uspp_F_UPF_bz2}
});

std::unordered_map<std::string, unsigned int> USPP_FILES_LEN ({
{"AG" ,  183881},
{"AL" ,  76968},
{"AS" ,  148775},
{"AU" ,  225819},
{"BA" ,  210009},
{"BE" ,  133580},
{"B" ,  85073},
{"BR" ,  149334},
{"CA" ,  152142},
{"CD" ,  175795},
{"CL" ,  152017},
{"CO" ,  158103},
{"C" ,  76929},
{"CR" ,  152100},
{"CS" ,  173991},
{"CU" ,  159196},
{"FE" ,  157983},
{"F" ,  77802},
{"GA" ,  157154},
{"GE" ,  154290},
{"HF" ,  250915},
{"HG" ,  185080},
{"H" ,  32685},
{"IN" ,  174892},
{"I" ,  171991},
{"IR" ,  229953},
{"K" ,  170681},
{"LA" ,  257695},
{"LI" ,  104960},
{"MG" ,  191579},
{"MN" ,  156358},
{"MO" ,  167710},
{"NA" ,  125524},
{"NB" ,  167829},
{"NI" ,  158757},
{"N" ,  77694},
{"O" ,  99143},
{"OS" ,  278489},
{"PB" ,  180288},
{"PD" ,  181370},
{"P" ,  152308},
{"PT" ,  189439},
{"RB" ,  162046},
{"RE" ,  280773},
{"RH" ,  220495},
{"RU" ,  184186},
{"SB" ,  174947},
{"SC" ,  158789},
{"SE" ,  148512},
{"SI" ,  152821},
{"SN" ,  175716},
{"S" ,  153061},
{"SR" ,  165425},
{"TA" ,  212553},
{"TC" ,  184433},
{"TE" ,  136023},
{"TI" ,  150704},
{"TL" ,  180838},
{"V" ,  152061},
{"W" ,  251353},
{"Y" ,  168179},
{"ZN" ,  174215},
{"ZR" ,  167538}
});
namespace io = boost::iostreams;

std::string GetInternalPseudo(const char *symbol)
{

    std::string decompressed;
    std::string lsym(symbol);
    boost::to_upper(lsym);

    char *pptr = (char *)USPP_FILES[lsym];
    unsigned int len = USPP_FILES_LEN[lsym];
    std::vector<char> pp_file_v(pptr, pptr + len);

    io::filtering_ostream os;
    os.push(io::bzip2_decompressor());
    os.push(io::back_inserter(decompressed));

    io::write(os, pptr, (std::streamsize)pp_file_v.size());
    return decompressed;

    
}

#endif
