/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"

extern double *vh_old, *vxc_old;

void init_TDDFT()
{
    double time1;
    int MAT_TRANSFER = 0;


    time1 = my_crtc();

    /* initialize processor structure, decompose the processor into subdomains
       pct.pe_kpoint * ( pct.pe_x, pct.pe_y, pct.pe_z) or
       pct.pe_kpoint * ( pct.pe_column , pct.pe_row)
     */
    init_dimension(&MXLLDA, &MXLCOL);
    init_pe_on();


    /* allocate memory for matrixs  */
    allocate_matrix();

    /* Perform some necessary initializations no matter localized or not  
     */
    my_malloc_init( vxc_old, get_FP0_BASIS(), double );
    my_malloc_init( vh_old, get_FP0_BASIS(), double );

    init(vh, rho, rhocore, rhoc, states, states1, vnuc, vxc, vh_old, vxc_old);

    my_barrier();


}
