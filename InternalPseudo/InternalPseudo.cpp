#if INTERNAL_PP

#include <vector>
#include <string>
#include "const.h"
#include "params.h"
#include "main.h"
#include "InternalPseudo.h"



std::string GetInternalPseudo(const char *symbol)
{
   if(ct.internal_pseudo_type == NORM_CONSERVING_SG15)
       return GetInternalPseudo_sg15(symbol);

   if(ct.internal_pseudo_type == ULTRASOFT_GBRV)
       return GetInternalPseudo_uspp(symbol);

   if(ct.internal_pseudo_type == NORM_CONSERVING_ACCURACY)
       return GetInternalPseudo_ncpp_stringent(symbol);

   if(ct.internal_pseudo_type == NORM_CONSERVING_STANDARD)
       return GetInternalPseudo_ncpp_standard(symbol);
}

#endif
