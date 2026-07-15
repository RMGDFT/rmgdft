#if INTERNAL_PP

#include <vector>
#include <string>
#include "const.h"
#include "params.h"
#include "main.h"
#include "InternalPseudo.h"


std::unordered_map<std::string, std::string> default_pp ({
{"Ag"  , "pd_standard"},
{"Al"  , "pd_standard"},
{"Ar"  , "pd_standard"},
{"As"  , "pd_standard"},
{"Au"  , "pd_standard"},
{"Ba"  , "pd_standard"},
{"Be"  , "pd_standard"},
{"Bi"  , "pd_standard"},
{"Br"  , "pd_standard"},
{"B"  , "pd_standard"},
{"Ca"  , "pd_standard"},
{"Cd"  , "pd_standard"},
{"Cl"  , "pd_standard"},
{"Co"  , "pd_standard"},
{"Cr"  , "pd_standard"},
{"Cs"  , "pd_standard"},
{"C"  , "pd_standard"},
{"Cu"  , "pd_standard"},
{"Fe"  , "pd_standard"},
{"F"  , "pd_standard"},
{"Ga"  , "pd_standard"},
{"Ge"  , "pd_standard"},
{"He"  , "pd_standard"},
{"Hf"  , "pd_standard"},
{"Hg"  , "pd_standard"},
{"H"  , "pd_standard"},
{"In"  , "pd_standard"},
{"Ir"  , "pd_standard"},
{"I"  , "pd_standard"},
{"Kr"  , "pd_standard"},
{"K"  , "pd_standard"},
{"La"  , "pd_standard"},
{"Li"  , "pd_standard"},
{"Lu"  , "pd_standard"},
{"Mg"  , "pd_standard"},
{"Mn"  , "pd_standard"},
{"Mo"  , "pd_standard"},
{"Na"  , "pd_standard"},
{"Nb"  , "pd_standard"},
{"Ne"  , "pd_standard"},
{"Ni"  , "pd_standard"},
{"N"  , "pd_standard"},
{"Os"  , "pd_standard"},
{"O"  , "pd_standard"},
{"Pb"  , "pd_standard"},
{"Pd"  , "pd_standard"},
{"Po"  , "pd_standard"},
{"P"  , "pd_standard"},
{"Pt"  , "pd_standard"},
{"Rb"  , "pd_standard"},
{"Re"  , "pd_standard"},
{"Rh"  , "pd_standard"},
{"Rn"  , "pd_standard"},
{"Ru"  , "pd_standard"},
{"Sb"  , "pd_standard"},
{"Sc"  , "pd_standard"},
{"Se"  , "pd_standard"},
{"Si"  , "pd_standard"},
{"Sn"  , "pd_standard"},
{"Sr"  , "pd_standard"},
{"S"  , "pd_standard"},
{"Ta"  , "pd_standard"},
{"Tc"  , "pd_standard"},
{"Te"  , "pd_standard"},
{"Ti"  , "pd_standard"},
{"Tl"  , "pd_standard"},
{"V"  , "pd_standard"},
{"W"  , "pd_standard"},
{"Xe"  , "pd_standard"},
{"Y"  , "pd_standard"},
{"Zn"  , "pd_standard"},
{"Zr"  , "pd_standard"}
});


std::string GetInternalPseudo(const char *symbol)
{
   if(ct.internal_pseudo_type == NORM_CONSERVING_DEFAULT)
   {
       std::string type = default_pp[std::string(symbol)];
       printf("RRRR  %s\n",type.c_str());
       if(type == std::string("sg15"))
           return GetInternalPseudo_sg15(symbol);
       if(type == std::string("nc_accuracy") || type == std::string("pd_precision"))
           return GetInternalPseudo_ncpp_stringent(symbol);
       if(type == std::string("nc_standard") || type == std::string("pd_standard"))
           return GetInternalPseudo_ncpp_standard(symbol);
   }

   if(ct.internal_pseudo_type == NORM_CONSERVING_SG15)
       return GetInternalPseudo_sg15(symbol);

   if(ct.internal_pseudo_type == ULTRASOFT_GBRV)
       return GetInternalPseudo_uspp(symbol);

   if(ct.internal_pseudo_type == NORM_CONSERVING_ACCURACY)
       return GetInternalPseudo_ncpp_stringent(symbol);

   if(ct.internal_pseudo_type == NORM_CONSERVING_STANDARD)
       return GetInternalPseudo_ncpp_standard(symbol);

   return std::string("");
}

#endif
