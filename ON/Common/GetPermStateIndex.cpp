
#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "rmgtypedefs.h"
#include "transition.h"

void GetPermStateIndex(int num_ions, ION *ions, unsigned int *perm_ion_index, 
unsigned int *perm_state_index, unsigned int *rev_perm_state_index)
{

    unsigned int *rev_perm;
    int ion, st_index, st, st_start;
    rev_perm = new unsigned int[num_ions];

    for(ion = 0; ion < num_ions; ion++)
        rev_perm[perm_ion_index[ion]] = ion;

    st_index = 0;
    for(ion = 0; ion < num_ions; ion++)
    {
        int ion_permuted = rev_perm[ion];
        st_start = ions[ion_permuted].orbital_start_index;

        for(st = 0; st < ions[ion_permuted].num_orbitals;st++)
            rev_perm_state_index[st_index + st] = st_start + st;
        st_index += ions[ion_permuted].num_orbitals;

    }

    for(st = 0; st < st_index; st++)
        perm_state_index[rev_perm_state_index[st]] = st;
        
    delete [] rev_perm;
} 

