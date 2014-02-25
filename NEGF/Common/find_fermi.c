/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


*/




#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "method.h"
#include "pmo.h"



void find_fermi (complex double * sigma_all)
{

    int fermi_step, st1, st2;
    rmg_double_t tem1, tem2, bias1, bias2, bias1a, bias2a, tchargea;
    rmg_double_t *matrix_SxRHO;
    int ntot, i;
    int ione = 1;

    ntot = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
    }

    for (i = 1; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i - 1] * pmo.mxlocc_cond[i];
    }

    my_malloc_init( matrix_SxRHO, ntot, rmg_double_t );

    for (fermi_step = 0; fermi_step < 4; fermi_step++)
    {
        if (fermi_step == 0)
        {
            bias1 = lcr[1].EF_new;
            bias2 = lcr[2].EF_new;
        }

        set_energy_weight (lcr[1].ene, lcr[1].weight, bias1, &lcr[1].nenergy);
        sigma_all_energy_point (sigma_all);

        charge_density_matrix_p (sigma_all);

        for (i = 0; i < ntot; i++)
            matrix_SxRHO[i] = lcr[0].density_matrix_tri[i] * lcr[0].Stri[i];

        /* calculating the total charge in the central part  */
        ct.tcharge = 0.0;
        for (st1 = pmo.offdiag_begin[0]; st1 < pmo.diag_begin[ct.num_blocks-1]; st1++)
                ct.tcharge += 6.0 * matrix_SxRHO[st1];
         comm_sums(&ct.tcharge, &ione, COMM_EN2);


        if (pct.gridpe == 0)
            printf ("\n total charge, %f %f \n Fermi energy %16.10f %16.10f",
                    ct.nel, ct.tcharge, bias1, bias2);


        if (fermi_step == 0)
        {
            bias1a = lcr[1].EF_new;
            bias2a = lcr[2].EF_new;
            tchargea = ct.tcharge;
            if (ct.nel > ct.tcharge)
            {
                bias1 += 0.05;
                bias2 += 0.05;
            }
            else
            {
                bias1 -= 0.05;
                bias2 -= 0.05;
            }
        }
        else
        {
            tem1 = bias1a - (bias1a - bias1) / (tchargea - ct.tcharge) * (tchargea - ct.nel);
            tem2 = bias2a - (bias2a - bias2) / (tchargea - ct.tcharge) * (tchargea - ct.nel);
            tchargea = ct.tcharge;
            bias1a = bias1;
            bias2a = bias2;
            bias1 = tem1;
            bias2 = tem2;
        }


    }

    if (pct.gridpe == 0)
        printf("\nFERMI ENERGY = %15.8f\n", bias1);

    my_free(matrix_SxRHO);
    set_energy_weight (lcr[1].ene, lcr[1].weight, bias1, &lcr[1].nenergy);

    sigma_all_energy_point (sigma_all);

    lcr[1].EF_new = bias1;
    lcr[2].EF_new = bias2;
    lcr[1].EF_old = bias1;
    lcr[2].EF_old = bias2;

}
