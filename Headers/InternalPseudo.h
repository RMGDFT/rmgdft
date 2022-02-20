#ifndef RMG_InternalPseudo_H
#define RMG_InternalPseudo_H 1

typedef unsigned char * compressed_pp;

std::string GetInternalPseudo(const char *symbol);
std::string GetInternalPseudo_uspp(const char *symbol);
std::string GetInternalPseudo_ncpp_stringent(const char *symbol);

#endif
