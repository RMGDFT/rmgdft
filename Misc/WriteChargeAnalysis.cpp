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
#include "RmgSumAll.h"
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



    rmg_printf ("\n\n");

    /*Get atomic charge density*/
    switch (ct.charge_analysis_type)
    {
	case CHARGE_ANALYSIS_NONE:
	    rmg_printf("No charge analysis performed (should not happen)");
	    break;

	case CHARGE_ANALYSIS_VORONOI:
	    rmg_printf("VORONOI DEFORMATION DENSITY");
	    break;

	default :
	    printf("Invalid Charge Analysis" );
    }


    rmg_printf (" PARTIAL CHARGES \n\n");
    rmg_printf("      Ion  Species      Charge\n");

    for (i = 0; i < ct.num_ions; i++)
    {

	iptr = &ct.ions[i];
	sp = &ct.sp[iptr->species];

	rmg_printf("     %3d     %2s          %6.3f\n", i + 1, sp->atomic_symbol, iptr->partial_charge); 

    }





}
