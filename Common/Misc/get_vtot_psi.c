/************************** SVN Revision Information **************************
 **    $Id: get_vtot_psi.c 1166 2010-11-17 20:22:00Z btan $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void get_vtot_psi (REAL * vtot_psi, REAL * vtot, int scale)
{

    int idx, ione =1;
    REAL *sg_vtot;

    if(scale ==1)
    {
        idx = FPX0_GRID * FPY0_GRID * FPZ0_GRID;    
        dcopy(&idx, vtot, &ione, vtot_psi, &ione);
        return;
    }

    my_malloc (sg_vtot,(FPX0_GRID+10)*(FPY0_GRID+10)*(FPZ0_GRID+10),REAL);
    trade_imagesx (vtot,sg_vtot,FPX0_GRID,FPY0_GRID,FPZ0_GRID,5);

    if(scale !=1 ) 
    {
        mg_restrict_6 (sg_vtot,vtot_psi,FPX0_GRID,FPY0_GRID,FPZ0_GRID,scale);
    } 
        

    my_free (sg_vtot);

}
