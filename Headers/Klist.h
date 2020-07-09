#ifndef KLIST_H
#define KLIST_H 1
#include <vector>
#include <complex>
#include <boost/multi_array.hpp>
// all kpoint infomations, whole kpoint, irreducible kpoints and thier symmetry relation

class Klist
{
public:
    int kpoint_mesh[3];
    int kpoint_is_shift[3];
    int num_k_all;
    int num_k_ire;
    
    boost::multi_array<double, 2> k_all_xtal;
    boost::multi_array<double, 2> k_ire_xtal;
    boost::multi_array<double, 2> k_all_cart;
    boost::multi_array<double, 2> k_ire_cart;

    std::vector<double> kweight;

    // all kpoint map to ireeducible kpoints.
    std::vector<int> k_map_index;
    std::vector<int> k_map_symm;

    boost::multi_array<int, 2> k_neighbor;
};
#endif
