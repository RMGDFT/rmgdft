#ifndef RMG_MapElements_H
#define RMG_MapElements_H 1
#include <string>

extern "C" int GetAtomicNumber(std::string symbol);
extern "C" const char * GetAtomicSymbol(int number);
extern "C" double GetAtomicMass(std::string symbol);

#endif
