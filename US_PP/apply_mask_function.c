/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <math.h>

/*Applies mask function
 *   f: a radial function to be filtered
 *   r: array stored r values
 *   rg_points: number of points in radial grid
 *   rcut: The distance after which f should be 0 (after filtering)
 *   offset: distance from the end at which the result is set to 0
 *           This function divides the original function by a mask function, which is zero at rcut. 
 *           For mask function to work correctly, f should be zero at rcut, but sometimes (especially for beta functions)
 *           f decays very slowly only to reach 0 at infinity. In order to avoid problems with dividing by zero we assume
 *           that at the distance (rcut-offset) f is zero and set f/mask to zero for r> (rcut-offset).
 * */           

#define VERBOSE 0


void apply_mask_function (double *f, double * r, int rg_points, double rmax, double offset)
{

int idx;
#if VERBOSE
    printf("\n\n");
#endif

    for (idx = 0; idx < rg_points; idx++)
    {   
	if (r[idx] < rmax-offset)
	    f[idx] = f[idx]/mask_function(r[idx]/rmax);
	else
            f[idx] = 0.0; 

#if VERBOSE
        printf("\n %f %f", r[idx], f[idx]);
#endif
    
}   


}

