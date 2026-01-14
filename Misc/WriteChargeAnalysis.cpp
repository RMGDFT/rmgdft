#include <vector>
#include <cstdlib>
#include <cmath>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "transition.h"
#include "rmg_sum_all.h"
#include "GlobalSums.h"
#include "boundary_conditions.h"


/*we can consider all ions (expensive) or just those nearby (may not be always sufficient)
 * Using local ions only, error is expected to be very small */

using namespace std;

void WriteChargeAnalysis(void)
{

    int i;
    ION *iptr;
    SPECIES *sp;



    rmg::printlog ("\n\n");

    /*Get atomic charge density*/
    switch (ct.charge_analysis_type)
    {
	case CHARGE_ANALYSIS_NONE:
	    rmg::printlog("No charge analysis performed (should not happen)");
	    break;

	case CHARGE_ANALYSIS_VORONOI:
	    rmg::printlog("VORONOI DEFORMATION DENSITY");
	    break;

	default :
	    printf("Invalid Charge Analysis" );
    }


    rmg::printlog (" PARTIAL CHARGES \n\n");
    rmg::printlog("      Ion  Species      Charge\n");

    for (i = 0; i < ct.num_ions; i++)
    {

	iptr = &Atoms[i];
	sp = &Species[iptr->species];

	rmg::printlog("     %3d     %2s          %6.3f\n", i + 1, sp->atomic_symbol, iptr->partial_charge); 

    }





}
