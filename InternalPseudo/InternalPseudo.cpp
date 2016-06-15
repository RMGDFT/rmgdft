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
namespace io = boost::iostreams;

std::string GetInternalPseudo(const char *symbol)
{

    std::string decompressed;
    std::string lsym(symbol);
    boost::to_upper(lsym);

    std::vector<unsigned char> pp_file_v = USPP_FILES[lsym];
    char *pptr = (char *)&pp_file_v[0];

    io::filtering_ostream os;
    os.push(io::bzip2_decompressor());
    os.push(io::back_inserter(decompressed));

    io::write(os, pptr, (std::streamsize)pp_file_v.size());
    return decompressed;

    
}

