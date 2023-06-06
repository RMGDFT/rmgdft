#ifndef RMG_MapElements_H
#define RMG_MapElements_H 1
#include <string>

extern "C" int GetAtomicNumber(std::string symbol);
extern "C" const char * GetAtomicSymbol(int number);
extern "C" double GetAtomicMass(std::string symbol);
int GetNumberOrbitalsL(std::string symbol);
int GetNumberOrbitalsM(std::string symbol);
void SetupAllElectonOrbitals(std::string symbol,
                             std::vector<int> &lvals,
                             std::vector<double> &jvals,
                             std::vector<double> &occs,
                             std::vector<double> &energy,
                             std::vector<double> &aradius,
                             std::vector<std::string> &label);


#endif
