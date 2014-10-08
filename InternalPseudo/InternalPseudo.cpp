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
{"H", h_pbe_v1_uspp_F_UPF_bz2},
{"C", c_pbe_v1_2_uspp_F_UPF_bz2},
{"N", n_pbe_v1_2_uspp_F_UPF_bz2},
{"O", o_pbe_v1_2_uspp_F_UPF_bz2}
});

namespace io = boost::iostreams;

std::string GetInternalPseudo(const char *symbol)
{

    std::string decompressed;
    std::string lsym(symbol);

    std::vector<unsigned char> pp_file_v = USPP_FILES[lsym];
    char *pptr = (char *)&pp_file_v[0];

    io::filtering_ostream os;
    os.push(io::bzip2_decompressor());
    os.push(io::back_inserter(decompressed));

    io::write(os, pptr, (std::streamsize)pp_file_v.size());
    return decompressed;

    
}

