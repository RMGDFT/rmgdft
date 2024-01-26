#include <unordered_map>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "MapElements.h"
#include "rmg_error.h"
#include "RmgException.h"
#include <iostream>


std::unordered_map<std::string, const int> SymbolToNumber ({
{"0?", 0},
{"H", 1},
{"He", 2},
{"Li", 3},
{"Be", 4},
{"B", 5},
{"C", 6},
{"N", 7},
{"O", 8},
{"F", 9},
{"Ne", 10},
{"Na", 11},
{"Mg", 12},
{"Al", 13},
{"Si", 14},
{"P", 15},
{"S", 16},
{"Cl", 17},
{"Ar", 18},
{"K", 19},
{"Ca", 20},
{"Sc", 21},
{"Ti", 22},
{"V", 23},
{"Cr", 24},
{"Mn", 25},
{"Fe", 26},
{"Co", 27},
{"Ni", 28},
{"Cu", 29},
{"Zn", 30},
{"Ga", 31},
{"Ge", 32},
{"As", 33},
{"Se", 34},
{"Br", 35},
{"Kr", 36},
{"Rb", 37},
{"Sr", 38},
{"Y", 39},
{"Zr", 40},
{"Nb", 41},
{"Mo", 42},
{"Tc", 43},
{"Ru", 44},
{"Rh", 45},
{"Pd", 46},
{"Ag", 47},
{"Cd", 48},
{"In", 49},
{"Sn", 50},
{"Sb", 51},
{"Te", 52},
{"I", 53},
{"Xe", 54},
{"Cs", 55},
{"Ba", 56},
{"La", 57},
{"Ce", 58},
{"Pr", 59},
{"Nd", 60},
{"Pm", 61},
{"Sm", 62},
{"Eu", 63},
{"Gd", 64},
{"Tb", 65},
{"Dy", 66},
{"Ho", 67},
{"Er", 68},
{"Tm", 69},
{"Yb", 70},
{"Lu", 71},
{"Hf", 72},
{"Ta", 73},
{"W", 74},
{"Re", 75},
{"Os", 76},
{"Ir", 77},
{"Pt", 78},
{"Au", 79},
{"Hg", 80},
{"Tl", 81},
{"Pb", 82},
{"Bi", 83},
{"Po", 84},
{"At", 85},
{"Rn", 86},
{"Fr", 87},
{"Ra", 88},
{"Ac", 89},
{"Th", 90},
{"Pa", 91},
{"U", 92},
{"Np", 93},
{"Pu", 94},
{"Am", 95},
{"Cm", 96},
{"Bk", 97},
{"Cf", 98},
{"Es", 99},
{"Fm", 100},
{"Md", 101},
{"No", 102},
{"Lr", 103},
{"Rf", 104},
{"Db", 105},
{"Sg", 106},
{"Bh", 107},
{"Hs", 108},
{"Mt", 109},
{"Ds", 110}
});

std::unordered_map<std::string, const int> SymbolToMass ({
{"H",  1.0080},
{"He",  4.0026},
{"Li",  6.9400},
{"Be",  9.0122},
{"B", 10.8100},
{"C", 12.0110},
{"N", 14.0070},
{"O", 15.9990},
{"F", 18.9984},
{"Ne", 20.1797},
{"Na", 22.9898},
{"Mg", 24.3050},
{"Al", 26.9815},
{"Si", 28.0850},
{"P", 30.9738},
{"S", 32.0600},
{"Cl", 35.4500},
{"Ar", 39.9480},
{"K", 39.0983},
{"Ca", 40.0780},
{"Sc", 44.9559},
{"Ti", 47.8670},
{"V", 50.9415},
{"Cr", 51.9961},
{"Mn", 54.9380},
{"Fe", 55.8450},
{"Co", 58.9332},
{"Ni", 58.6934},
{"Cu", 63.5460},
{"Zn", 65.3800},
{"Ga", 69.7230},
{"Ge", 72.6300},
{"As", 74.9216},
{"Se", 78.9710},
{"Br", 79.9040},
{"Kr", 83.7980},
{"Rb", 85.4678},
{"Sr", 87.6200},
{"Y", 88.9058},
{"Zr", 91.2240},
{"Nb", 92.9064},
{"Mo", 95.9500},
{"Tc", 97.0000},
{"Ru", 101.0700},
{"Rh", 102.9055},
{"Pd", 106.4200},
{"Ag", 107.8682},
{"Cd", 112.4140},
{"In", 114.8180},
{"Sn", 118.7100},
{"Sb", 121.7600},
{"Te", 127.6000},
{"I", 126.9045},
{"Xe", 131.2930},
{"Cs", 132.9055},
{"Ba", 137.3270},
{"La", 138.9055},
{"Ce", 140.1160},
{"Pr", 140.9077},
{"Nd", 144.2420},
{"Pm", 145.0000},
{"Sm", 150.3600},
{"Eu", 151.9640},
{"Gd", 157.2500},
{"Tb", 158.9254},
{"Dy", 162.5000},
{"Ho", 164.9303},
{"Er", 167.2590},
{"Tm", 168.9342},
{"Yb", 173.0540},
{"Lu", 174.9668},
{"Hf", 178.4900},
{"Ta", 180.9479},
{"W", 183.8400},
{"Re", 186.2070},
{"Os", 190.2300},
{"Ir", 192.2170},
{"Pt", 195.0840},
{"Au", 196.9666},
{"Hg", 200.5920},
{"Tl", 204.3800},
{"Pb", 207.2000},
{"Bi", 208.9804},
{"Po", 209.0000},
{"At", 210.0000},
{"Rn", 222.0000},
{"Fr", 223.0000},
{"Ra", 226.0000},
{"Ac", 227.0000},
{"Th", 232.0377},
{"Pa", 231.0359},
{"U", 238.0289},
{"Np", 237.0000},
{"Pu", 244.0000},
{"Am", 243.0000},
{"Cm", 247.0000},
{"Bk", 247.0000},
{"Cf", 251.0000},
{"Es", 252.0000},
{"Fm", 257.0000},
{"Md", 258.0000},
{"No", 259.0000},
{"Lr", 262.0000},
{"Rf", 267.0000},
{"Db", 270.0000},
{"Sg", 271.0000},
{"Bh", 270.0000},
{"Hs", 277.0000},
{"Mt", 276.0000},
{"Ds", 281.0000},
{"DLO", 1000000.0000},
});


std::unordered_map<std::string, std::string> SymbolToConfig ({
{"H", "1s1"},
{"He", "1s2"},
{"Li", "1s2 2s1"},
{"Be", "1s2 2s2"},
{"B", "1s2 2s2 2p1"},
{"C", "1s2 2s2 2p2"},
{"N", "1s2 2s2 2p3"},
{"O", "1s2 2s2 2p4"},
{"F", "1s2 2s2 2p5"},
{"Ne", "1s2 2s2 2p6"},
{"Na", "1s2 2s2 2p6 3s1"},
{"Mg", "1s2 2s2 2p6 3s2"},
{"Al", "1s2 2s2 2p6 3s2 3p1"},
{"Si", "1s2 2s2 2p6 3s2 3p2"},
{"P", "1s2 2s2 2p6 3s2 3p3"},
{"S", "1s2 2s2 2p6 3s2 3p4"},
{"Cl", "1s2 2s2 2p6 3s2 3p5"},
{"Ar", "1s2 2s2 2p6 3s2 3p6"},
{"K", "1s2 2s2 2p6 3s2 3p6 4s1"},
{"Ca", "1s2 2s2 2p6 3s2 3p6 4s2"},
{"Sc", "1s2 2s2 2p6 3s2 3p6 3d1 4s2"},
{"Ti", "1s2 2s2 2p6 3s2 3p6 3d2 4s2"},
{"V", "1s2 2s2 2p6 3s2 3p6 3d3 4s2"},
{"Cr", "1s2 2s2 2p6 3s2 3p6 3d5 4s1"},
{"Mn", "1s2 2s2 2p6 3s2 3p6 3d5 4s2"},
{"Fe", "1s2 2s2 2p6 3s2 3p6 3d6 4s2"},
{"Co", "1s2 2s2 2p6 3s2 3p6 3d7 4s2"},
{"Ni", "1s2 2s2 2p6 3s2 3p6 3d8 4s2"},
{"Cu", "1s2 2s2 2p6 3s2 3p6 3d10 4s1"},
{"Zn", "1s2 2s2 2p6 3s2 3p6 3d10 4s2"},
{"Ga", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p1"},
{"Ge", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p2"},
{"As", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p3"},
{"Se", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p4"},
{"Br", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p5"},
{"Kr", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6"},
{"Rb", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s1"},
{"Sr", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2"},
{"Y", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d1 5s2"},
{"Zr", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d2 5s2"},
{"Nb", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d4 5s1"},
{"Mo", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d5 5s1"},
{"Tc", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d5 5s2"},
{"Ru", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d7 5s1"},
{"Rh", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d8 5s1"},
{"Pd", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10"},
{"Ag", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s1"},
{"Cd", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2"},
{"In", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p1"},
{"Sn", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p2"},
{"Sb", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p3"},
{"Te", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p4"},
{"I", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p5"},
{"Xe", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6"},
{"Cs", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 6s1"},
{"Ba", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 6s2"},
{"La", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 5d1 6s2"},
{"Ce", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f1 5s2 5p6 5d1 6s2"},
{"Pr", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f3 5s2 5p6 6s2"},
{"Nd", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f4 5s2 5p6 6s2"},
{"Pm", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f5 5s2 5p6 6s2"},
{"Sm", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f6 5s2 5p6 6s2"},
{"Eu", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f7 5s2 5p6 6s2"},
{"Gd", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f7 5s2 5p6 5d1 6s2"},
{"Tb", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f9 5s2 5p6 6s2"},
{"Dy", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f10 5s2 5p6 6s2"},
{"Ho", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f11 5s2 5p6 6s2"},
{"Er", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f12 5s2 5p6 6s2"},
{"Tm", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f13 5s2 5p6 6s2"},
{"Yb", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 6s2"},
{"Lu", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d1 6s2"},
{"Hf", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d2 6s2"},
{"Ta", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d3 6s2"},
{"W", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d4 6s2"},
{"Re", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d5 6s2"},
{"Os", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d6 6s2"},
{"Ir", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d7 6s2"},
{"Pt", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d9 6s1"},
{"Au", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s1"},
{"Hg", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2"},
{"Tl", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p1"},
{"Pb", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p2"},
{"Bi", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p3"},
{"Po", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p4"},
{"At", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p5"},
{"Rn", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6"},
{"Fr", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6 7s1"},
{"Ra", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6 7s2"},
{"Ac", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6 6d1 7s2"},
{"Th", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6 6d2 7s2"},
{"Pa", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f2 6s2 6p6 6d1 7s2"},
{"U", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f3 6s2 6p6 6d1 7s2"},
{"Np", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f4 6s2 6p6 6d1 7s2"},
{"Pu", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f6 6s2 6p6 7s2"},
{"Am", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f7 6s2 6p6 7s2"},
{"Cm", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f7 6s2 6p6 6d1 7s2"},
{"Bk", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f9 6s2 6p6 7s2"},
{"Cf", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f10 6s2 6p6 7s2"},
{"Es", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f11 6s2 6p6 7s2"},
{"Fm", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f12 6s2 6p6 7s2"},
{"Md", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f13 6s2 6p6 7s2"},
{"No", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 7s2"},
{"Lr", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 7s2 7p1"},
{"Rf", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 6d2 7s2"},
{"Db", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 6d3 7s2"},
{"Sg", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 6d4 7s2"},
{"Bh", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 6d5 7s2"},
{"Hs", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 6d6 7s2"},
{"Mt", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 6d7 7s2"},
{"Ds", "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 5f14 6s2 6p6 6d8 7s2"},
});

void split(std::string const &str, const char token,
            std::vector<std::string> &out)
{
    size_t begin, end=0;
 
    while ((begin = str.find_first_not_of(token, end)) != std::string::npos)
    {
        end = str.find(token, begin);
        out.push_back(str.substr(begin, end - begin));
    }
}

void CheckSymbol(std::string const &symbol)
{
    std::string ss;
    try {
        ss = SymbolToConfig[symbol];
    }
    catch (const std::out_of_range& oor) {
        throw RmgFatalException() << "Unknown atomic symbol " << symbol << " in " << __FILE__ << " at line " << __LINE__ << "\n";
    }
}

int GetAtomicNumber(std::string symbol)
{
    CheckSymbol(symbol);
    return SymbolToNumber[symbol];
}

// Counts the number of radial atomic orbitals 
int GetNumberOrbitalsL(std::string symbol)
{
    CheckSymbol(symbol);
    std::string ss;
    ss = SymbolToConfig[symbol];
    const char token = ' ';
    std::vector<std::string> fields;
    split(ss, token, fields);
    return (int)fields.size();
}

// Counts the number of atomic orbitals including angular dependence
int GetNumberOrbitalsM(std::string symbol)
{
    CheckSymbol(symbol);
    std::string ss;
    ss = SymbolToConfig[symbol];

    int norbs = 0;
    const char token = ' ';
    std::vector<std::string> fields;
    split(ss, token, fields);
    for (auto &field: fields)
    {
        if(field.find('s') != std::string::npos) norbs = norbs + 1;
        if(field.find('p') != std::string::npos) norbs = norbs + 3;
        if(field.find('d') != std::string::npos) norbs = norbs + 5;
        if(field.find('f') != std::string::npos) norbs = norbs + 7;
    }
    return norbs;
}

double GetAtomicMass(std::string symbol)
{
    CheckSymbol(symbol);
    return SymbolToMass[symbol];
}

const char * GetAtomicSymbol(int number)
{
    if((number < 1) || (number > 110))
        throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << ". Atomic numbers must be between 1 and 110.\n";

    for ( auto it = SymbolToNumber.begin(); it != SymbolToNumber.end(); ++it ) {
        if(it->second == number) {
            return it->first.c_str();
        }
    }
    return "";
}    

void SetupAllElectonOrbitals(std::string symbol, 
                             std::vector<int> &pqn, 
                             std::vector<int> &lvals, 
                             std::vector<double> &jvals,
                             std::vector<double> &occs,
                             std::vector<double> &energy,
                             std::vector<double> &aradius,
                             std::vector<std::string> &label)
{
    CheckSymbol(symbol);
    std::string ss;
    ss = SymbolToConfig[symbol];

    const char token = ' ';
    std::vector<std::string> fields;
    split(ss, token, fields);
    for (auto &field: fields)
    {
        label.emplace_back(field);
        int p1 = stoi(field.substr(0,1));
        pqn.emplace_back(p1);
        double occ = stod(field.substr(2));
        occs.emplace_back(occ);
        energy.emplace_back(0.0);
        aradius.emplace_back(12.0);
        jvals.emplace_back(0);  // Nothing here now
        if(field.find('s') != std::string::npos)
        {
            lvals.emplace_back(0);
        }
        else if(field.find('p') != std::string::npos)
        {
            lvals.emplace_back(1);
        }
        else if(field.find('d') != std::string::npos)
        {
            lvals.emplace_back(2);
        }
        else if(field.find('f') != std::string::npos)
        {
            lvals.emplace_back(3);
        }
        else if(field.find('g') != std::string::npos)
        {
            lvals.emplace_back(4);
        }
    }
}
