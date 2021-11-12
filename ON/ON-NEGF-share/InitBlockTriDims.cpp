#include <float.h>
#include <stdio.h>
#include <assert.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"
#include "LocalObject.h"
#include "blacs.h"
#include "RmgException.h"
void InitBlockTriDims()
{
    double max_orbital_radius = 0.0;
    for(int st = 0; st < ct.num_states; st++) max_orbital_radius = std::max(max_orbital_radius, states[st].radius);
    double max_nl_radius = 0.0;
    for(int sp = 0; sp < ct.num_species; sp++) max_nl_radius = std::max(max_nl_radius, Species[sp].nlradius);
    double nbl =Rmg_L.celldm[0] /(max_orbital_radius + max_nl_radius)/2.0;

    ct.num_blocks = (int)nbl;
    //printf("\n max radius %f %f %f %d\n", max_orbital_radius, max_nl_radius, nbl, ct.num_blocks );
    for(int st = 1; st < ct.num_states; st++) {
        if( (states[st].crds[0] - states[st-1].crds[0]) <-1.0e-5)
           throw RmgFatalException() << "orbital x coordinates not sorted  \n";
 
    }

    ct.block_dim_phi.resize(ct.num_blocks, 0);
    ct.block_dim_nl.resize(ct.num_blocks, 0);

    double block_length_x = Rmg_L.celldm[0]/ct.num_blocks;
    for(int st = 0; st < ct.num_states; st++) {
        int ib = (int)(states[st].crds[0]/block_length_x );
        ct.block_dim_phi[ib]++;
    }

    for(int ion = 0; ion < ct.num_ions; ion++) {
        ION *iptr = &Atoms[ion];
        SPECIES *sp = &Species[iptr->species];
        int nh = sp->num_projectors;

        int ib = (int)(Atoms[ion].crds[0]/block_length_x );
        ct.block_dim_nl[ib] += nh;
        
    }

    rmg_printf("\n number of blocks %d and block dims for orbital and non-local projector", ct.num_blocks);
    for(int ib = 0; ib <ct.num_blocks; ib++) rmg_printf("\n block %d %d %d", ib, ct.block_dim_phi[ib], ct.block_dim_nl[ib]);
}
