/************************** SVN Revision Information **************************
 **    $Id: 
******************************************************************************/

/* calculating Hamiltonian matrix Hij and 
 * overlap matrix matB togather
 */

 
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include "make_conf.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "blas.h"
#include "Kbpsi.h"



void KbpsiUpdate(STATE * states)
{
    int idx, st1;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    int ixx, iyy, izz;
    unsigned int ion, num_orbital_thision, num_proj;
    int ip, iip1;

    RmgTimer *RT2 = new RmgTimer("2-SCF: kbpsi: calc");
    get_all_kbpsi(states, states, ion_orbit_overlap_region_nl, projectors, kbpsi);
    delete(RT2);

    for(ion = 0; ion < pct.n_ion_center; ion++)
    {
        Kbpsi_str.kbpsi_ion[ion].clear();
        num_orbital_thision = Kbpsi_str.num_orbital_thision[ion];
        num_proj = pct.prj_per_ion[pct.ionidx[ion]];
        Kbpsi_str.orbital_index[ion].resize(num_orbital_thision);


        // uopdate values from this process
        for(idx = 0; idx < num_orbital_thision; idx++)
        {
            st1 =Kbpsi_str.orbital_index[ion][idx];
            iip1 = (st1-ct.state_begin) * pct.n_ion_center * ct.max_nl + ion * ct.max_nl;
            for(ip = 0; ip < num_proj; ip++)
                Kbpsi_str.kbpsi_ion[ion].emplace_back(kbpsi[iip1 + ip]);
        }


    }

    RmgTimer *RT2a = new RmgTimer("2-SCF: kbpsi: comm");
    KbpsiComm();
    delete(RT2a);

}
