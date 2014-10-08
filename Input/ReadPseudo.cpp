#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "MapElements.h"
#include "transition.h"


extern "C" void ReadPseudo(int nspecies, SPECIES *sp)
{

    // Default to norm conserving pseudopotentials. If any of the atomic species
    // are represented by ultrasoft then all of them must be handled as ultrasoft
    ct.norm_conserving_pp = true;

    for(int isp = 0;isp < nspecies;isp++) {
        LoadUpf(sp);
        if(!sp->is_norm_conserving) ct.norm_conserving_pp = false;
        sp++; 
    } 


    // Now verify the pseudopotentials. Specifically all pseudopotentials must use the same
    // XC type. If any disagree then the value must be set in the input file and will override
    // the settings in the individual pseudopotentials.
}
