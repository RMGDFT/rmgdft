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
#include "pseudo_list_sg15.h"

std::unordered_map<std::string, compressed_pp> SG15_FILES ({
{"AG" ,	Ag_sg15_upf_bz2},
{"AL" ,	Al_sg15_upf_bz2},
{"AR" ,	Ar_sg15_upf_bz2},
{"AS" ,	As_sg15_upf_bz2},
{"AU" ,	Au_sg15_upf_bz2},
{"BA" ,	Ba_sg15_upf_bz2},
{"BE" ,	Be_sg15_upf_bz2},
{"BI" ,	Bi_sg15_upf_bz2},
{"BR" ,	Br_sg15_upf_bz2},
{"B" ,	B_sg15_upf_bz2},
{"CA" ,	Ca_sg15_upf_bz2},
{"CD" ,	Cd_sg15_upf_bz2},
{"CL" ,	Cl_sg15_upf_bz2},
{"CO" ,	Co_sg15_upf_bz2},
{"CR" ,	Cr_sg15_upf_bz2},
{"C" ,	C_sg15_upf_bz2},
{"CS" ,	Cs_sg15_upf_bz2},
{"CU" ,	Cu_sg15_upf_bz2},
{"FE" ,	Fe_sg15_upf_bz2},
{"F" ,	F_sg15_upf_bz2},
{"GA" ,	Ga_sg15_upf_bz2},
{"GE" ,	Ge_sg15_upf_bz2},
{"HE" ,	He_sg15_upf_bz2},
{"HF" ,	Hf_sg15_upf_bz2},
{"HG" ,	Hg_sg15_upf_bz2},
{"H" ,	H_sg15_upf_bz2},
{"IN" ,	In_sg15_upf_bz2},
{"IR" ,	Ir_sg15_upf_bz2},
{"I" ,	I_sg15_upf_bz2},
{"KR" ,	Kr_sg15_upf_bz2},
{"K" ,	K_sg15_upf_bz2},
{"LA" ,	La_sg15_upf_bz2},
{"LI" ,	Li_sg15_upf_bz2},
{"MG" ,	Mg_sg15_upf_bz2},
{"MN" ,	Mn_sg15_upf_bz2},
{"MO" ,	Mo_sg15_upf_bz2},
{"NA" ,	Na_sg15_upf_bz2},
{"NB" ,	Nb_sg15_upf_bz2},
{"NE" ,	Ne_sg15_upf_bz2},
{"NI" ,	Ni_sg15_upf_bz2},
{"N" ,	N_sg15_upf_bz2},
{"O" ,	O_sg15_upf_bz2},
{"OS" ,	Os_sg15_upf_bz2},
{"PB" ,	Pb_sg15_upf_bz2},
{"PD" ,	Pd_sg15_upf_bz2},
{"PO" ,	Po_sg15_upf_bz2},
{"P" ,	P_sg15_upf_bz2},
{"PT" ,	Pt_sg15_upf_bz2},
{"RB" ,	Rb_sg15_upf_bz2},
{"RE" ,	Re_sg15_upf_bz2},
{"RH" ,	Rh_sg15_upf_bz2},
{"RU" ,	Ru_sg15_upf_bz2},
{"SB" ,	Sb_sg15_upf_bz2},
{"SC" ,	Sc_sg15_upf_bz2},
{"SE" ,	Se_sg15_upf_bz2},
{"SI" ,	Si_sg15_upf_bz2},
{"SN" ,	Sn_sg15_upf_bz2},
{"SR" ,	Sr_sg15_upf_bz2},
{"S" ,	S_sg15_upf_bz2},
{"TA" ,	Ta_sg15_upf_bz2},
{"TC" ,	Tc_sg15_upf_bz2},
{"TE" ,	Te_sg15_upf_bz2},
{"TI" ,	Ti_sg15_upf_bz2},
{"TL" ,	Tl_sg15_upf_bz2},
{"V" ,	V_sg15_upf_bz2},
{"W" ,	W_sg15_upf_bz2},
{"XE" ,	Xe_sg15_upf_bz2},
{"Y" ,	Y_sg15_upf_bz2},
{"ZN" ,	Zn_sg15_upf_bz2},
{"ZR" ,	Zr_sg15_upf_bz2}
});

std::unordered_map<std::string, unsigned int> SG15_FILES_LEN ({
{"AG" , 120039},
{"AL" , 91882},
{"AR" , 80003},
{"AS" , 103029},
{"AU" , 143289},
{"BA" , 135206},
{"BE" , 73103},
{"BI" , 139938},
{"BR" , 103731},
{"B" , 73975},
{"CA" , 106057},
{"CD" , 119989},
{"CL" , 79974},
{"CO" , 114174},
{"CR" , 113468},
{"C" , 74025},
{"CS" , 115729},
{"CU" , 114496},
{"FE" , 114083},
{"F" , 75567},
{"GA" , 113613},
{"GE" , 112744},
{"HE" , 46948},
{"HF" , 157221},
{"HG" , 143595},
{"H" , 43937},
{"IN" , 116368},
{"IR" , 142649},
{"I" , 116559},
{"KR" , 102950},
{"K" , 106040},
{"LA" , 141599},
{"LI" , 69188},
{"MG" , 84090},
{"MN" , 114446},
{"MO" , 118822},
{"NA" , 83106},
{"NB" , 121040},
{"NE" , 77009},
{"NI" , 115639},
{"N" , 75651},
{"O" , 75286},
{"OS" , 144213},
{"PB" , 138755},
{"PD" , 119328},
{"PO" , 116612},
{"P" , 82075},
{"PT" , 142869},
{"RB" , 109677},
{"RE" , 143429},
{"RH" , 119506},
{"RU" , 119242},
{"SB" , 116275},
{"SC" , 112863},
{"SE" , 108193},
{"SI" , 81997},
{"SN" , 115180},
{"SR" , 109349},
{"S" , 81907},
{"TA" , 151761},
{"TC" , 118769},
{"TE" , 115738},
{"TI" , 113317},
{"TL" , 144408},
{"V" , 113862},
{"W" , 152232},
{"XE" , 115849},
{"Y" , 117747},
{"ZN" , 115746},
{"ZR" , 117296}
});
namespace io = boost::iostreams;

std::string GetInternalPseudo_sg15(const char *symbol)
{

    std::string decompressed;
    std::string lsym(symbol);
    boost::to_upper(lsym);

    char *pptr = (char *)SG15_FILES[lsym];
    unsigned int len = SG15_FILES_LEN[lsym];
    std::vector<char> pp_file_v(pptr, pptr + len);

    io::filtering_ostream os;
    os.push(io::bzip2_decompressor());
    os.push(io::back_inserter(decompressed));

    os.flush();
    io::write(os, pptr, (std::streamsize)pp_file_v.size());
    os.flush();
    io::close(os);
    // Little bit of a hack to remove some extraneous chars that were messing up the
    // UPF parsing. Fix later in pp generation.
//    decompressed.erase(decompressed.length() - 9);
    return decompressed;

    
}

#endif
