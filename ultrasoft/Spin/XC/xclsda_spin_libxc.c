
#include "main.h"
#include <float.h>
#include <math.h> 
#include "../../../Common/libxc/include/xc.h"   
/* include Libxc's header file */


void xclsda_spin_libxc (REAL * rho, REAL * rho_oppo, REAL * vxc, REAL * exc)
{

   /* XC_LDA_X = 1, XC_LDA_C_PZ = 9, XC_LDA_C_PW = 12*/
    int func_id_x = 1, func_id_c = 9; 
    xc_func_type func_x, func_c;
    /* double rhospin(*)[2], vspin(*)[2]; */
    REAL rhospin[FP0_BASIS][2], vspin[FP0_BASIS][2];
    REAL *ec;
    int idx; 

    	
    my_calloc (ec, FP0_BASIS, REAL);
    
    if(xc_func_init(&func_x, func_id_x, XC_POLARIZED) != 0)
    {
	if (pct.imgpe == 0)    
    		printf( "Functional %d not found\n", func_id_x);
    	exit(1);
    }
    
    if(xc_func_init(&func_c, func_id_c, XC_POLARIZED) != 0)
    {
	if (pct.imgpe == 0)    
    		printf( "Functional %d not found\n", func_id_c);
    	exit(1);
    } 


    /* pack rho and rho_oppo into the 2D array rhospin */ 
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    rhospin[idx][0] = rho[idx];
	    rhospin[idx][1] = rho_oppo[idx];
    }
    

    /* get the exchange part for energy and potential first*/    
    if (func_x.info->family == XC_FAMILY_LDA)
        xc_lda_exc_vxc (&func_x, FP0_BASIS, &rhospin[0][0], exc, &vspin[0][0]); 


    /* unpack the 2D array of exchange potential to vxc */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    vxc[idx] = vspin[idx][0];
    }
    

    /* get the correlation part for energy and potential */    
    if (func_c.info->family == XC_FAMILY_LDA)
        xc_lda_exc_vxc (&func_c, FP0_BASIS, &rhospin[0][0], ec, &vspin[0][0]); 


    /* add exchange and correlation together */
    for (idx = 0; idx < FP0_BASIS; idx++) 
    {
	    vxc[idx] += vspin[idx][0];
	    exc[idx] += ec[idx];
    }
    
    xc_func_end (&func_x);
    xc_func_end (&func_c);  


    my_free (ec);


}
