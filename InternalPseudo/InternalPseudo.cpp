#if INTERNAL_PP

#include <vector>
#include <string>
#include "const.h"
#include "params.h"
#include "main.h"
#include "InternalPseudo.h"


std::unordered_map<std::string, std::string> default_pp ({
{"AG"  , "abinit_standard"},
{"AL"  , "abinit_standard"},
{"AR"  , "abinit_standard"},
{"AS"  , "abinit_standard"},
{"AU"  , "abinit_standard"},
{"BA"  , "abinit_standard"},
{"BE"  , "abinit_standard"},
{"BI"  , "abinit_standard"},
{"BR"  , "abinit_standard"},
{"B"  , "abinit_standard"},
{"CA"  , "abinit_standard"},
{"CD"  , "abinit_standard"},
{"CL"  , "abinit_standard"},
{"CO"  , "abinit_standard"},
{"CR"  , "abinit_standard"},
{"CS"  , "abinit_standard"},
{"C"  , "abinit_standard"},
{"CU"  , "abinit_standard"},
{"FE"  , "abinit_standard"},
{"F"  , "abinit_standard"},
{"GA"  , "abinit_standard"},
{"GE"  , "abinit_standard"},
{"HE"  , "abinit_standard"},
{"HF"  , "abinit_standard"},
{"HG"  , "abinit_standard"},
{"H"  , "abinit_standard"},
{"IN"  , "abinit_standard"},
{"IR"  , "abinit_standard"},
{"I"  , "abinit_standard"},
{"KR"  , "abinit_standard"},
{"K"  , "abinit_standard"},
{"LA"  , "abinit_standard"},
{"LI"  , "abinit_standard"},
{"LU"  , "abinit_standard"},
{"MG"  , "abinit_standard"},
{"MN"  , "abinit_standard"},
{"MO"  , "abinit_standard"},
{"NA"  , "abinit_standard"},
{"NB"  , "abinit_standard"},
{"NE"  , "abinit_standard"},
{"NI"  , "abinit_standard"},
{"N"  , "abinit_standard"},
{"OS"  , "abinit_standard"},
{"O"  , "abinit_standard"},
{"PB"  , "abinit_standard"},
{"PD"  , "abinit_standard"},
{"PO"  , "abinit_standard"},
{"P"  , "abinit_standard"},
{"PT"  , "abinit_standard"},
{"RB"  , "abinit_standard"},
{"RE"  , "abinit_standard"},
{"RH"  , "abinit_standard"},
{"RN"  , "abinit_standard"},
{"RU"  , "abinit_standard"},
{"SB"  , "abinit_standard"},
{"SC"  , "abinit_standard"},
{"SE"  , "abinit_standard"},
{"SI"  , "abinit_standard"},
{"SN"  , "abinit_standard"},
{"SR"  , "abinit_standard"},
{"S"  , "abinit_standard"},
{"TA"  , "abinit_standard"},
{"TC"  , "abinit_standard"},
{"TE"  , "abinit_standard"},
{"TI"  , "abinit_standard"},
{"TL"  , "abinit_standard"},
{"V"  , "abinit_standard"},
{"W"  , "abinit_standard"},
{"XE"  , "abinit_standard"},
{"Y"  , "abinit_standard"},
{"ZN"  , "abinit_standard"},
{"ZR"  , "abinit_standard"}
});


std::string GetInternalPseudo(const char *symbol)
{
   if(ct.internal_pseudo_type == NORM_CONSERVING_DEFAULT)
   {
       std::string type = default_pp[symbol];
       if(type == std::string("sg15"))
           return GetInternalPseudo_sg15(symbol);
       if(type == std::string("nc_accuracy") || type == std::string("abinit_precision"))
           return GetInternalPseudo_ncpp_stringent(symbol);
       if(type == std::string("nc_standard") || type == std::string("abinit_standard"))
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
