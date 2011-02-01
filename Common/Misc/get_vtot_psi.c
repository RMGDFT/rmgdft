/************************** SVN Revision Information **************************
 **    $Id: get_vtot_psi.c 1166 2010-11-17 20:22:00Z btan $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void get_vtot_psi (REAL * vtot_psi, REAL * vtot)
{

   REAL *sg_vtot;

   my_malloc (sg_vtot,(FPX0_GRID+10)*(FPY0_GRID+10)*(FPZ0_GRID+10),REAL);
   trade_imagesx (vtot,sg_vtot,FPX0_GRID,FPY0_GRID,FPZ0_GRID,5);

   mg_restrict_6 (sg_vtot,vtot_psi,FPX0_GRID,FPY0_GRID,FPZ0_GRID);
   my_free (sg_vtot);

/*

    int i, j, k;       
    P0_GRID *ptr_vtot_psi;
    FP0_GRID *ptr_vtot;
 
    ptr_vtot_psi = (P0_GRID *) vtot_psi;
    ptr_vtot = (FP0_GRID *) vtot;


    for (i = 0; i < PX0_GRID; i++)
    {
        for (j = 0; j < PY0_GRID; j++)
        {
            for (k = 0; k < PZ0_GRID; k++)
            {
                ptr_vtot_psi->s1.b[i][j][k] = ptr_vtot->s1.b[FG_NX * i][FG_NY * j][FG_NZ * k];
            }
        }
    }
*/  


}
