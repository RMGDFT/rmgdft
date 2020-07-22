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
    int num_k_ext;
    int num_k_ire;
    int num_k_nn;
    
    boost::multi_array<double, 2> k_all_xtal;
    boost::multi_array<double, 2> k_ext_xtal;
    boost::multi_array<double, 2> k_ire_xtal;
    boost::multi_array<double, 2> k_all_cart;
    boost::multi_array<double, 2> k_ire_cart;

    std::vector<double> kweight;

    // all kpoint map to ireeducible kpoints.
    std::vector<int> k_map_index;
    std::vector<int> k_map_symm;

    //k_neightbors[num_k_all][num_k_nn][4];
    //[0]: kindex
    //[1,2,3], folding the kpoint into the real neighbor

    boost::multi_array<int, 3> k_neighbors;

};
#endif
