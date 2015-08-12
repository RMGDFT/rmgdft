
#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "rmgtypedefs.h"
#include "transition.h"

void PermAtoms(int num_ions, ION *ions, int *perm_index)
{

    ION *ions_perm;
    ions_perm = new ION[num_ions];
    for(int ion = 0; ion < num_ions; ion++)
    {
      ions_perm[ion] = ions[ion];
      //ions_perm[ion].crds[0] = ions[ion].crds[0];
      //ions_perm[ion].crds[1] = ions[ion].crds[1];
      //ions_perm[ion].crds[2] = ions[ion].crds[2];
      //ions_perm[ion].xtal[0] = ions[ion].xtal[0];
      //ions_perm[ion].xtal[1] = ions[ion].xtal[1];
      //ions_perm[ion].xtal[2] = ions[ion].xtal[2];
      //ions_perm[ion].species = ions[ion].species;
    }
    for(int ion = 0; ion < num_ions; ion++)
    {
      ions[ion] = ions_perm[perm_index[ion]];
      //ions[ion].crds[0] = ions_perm[perm_index[ion]].crds[0];
      //ions[ion].crds[1] = ions_perm[perm_index[ion]].crds[1];
      //ions[ion].crds[2] = ions_perm[perm_index[ion]].crds[2];
      //ions[ion].xtal[0] = ions_perm[perm_index[ion]].xtal[0];
      //ions[ion].xtal[1] = ions_perm[perm_index[ion]].xtal[1];
      //ions[ion].xtal[2] = ions_perm[perm_index[ion]].xtal[2];
      //ions[ion].species = ions_perm[perm_index[ion]].species;
    }
    delete [] ions_perm;
} 

