#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id: find_fermi.c 3140 2015-08-06 15:48:24Z luw $    **
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



void find_fermi (std::complex<double> * sigma_all)
{

    int fermi_step, st1;
    double bias1, bias1a=0.0, tchargea=0.0;
    int ione = 1;
    static double density = 0.0;


    bias1 = lcr[1].EF_new;
    for (fermi_step = 0; fermi_step < 4; fermi_step++)
    {

        set_energy_weight (lcr[1].ene, lcr[1].weight, bias1, &lcr[1].nenergy);
        sigma_all_energy_point (sigma_all, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2]);

        charge_density_matrix_p (sigma_all);

        /* calculating the total charge in the central part  */
        ct.tcharge = 0.0;
        for (st1 = pmo.offdiag_begin[0]; st1 < pmo.diag_begin[ct.num_blocks-1]; st1++)
                ct.tcharge += 6.0 * lcr[0].density_matrix_tri[st1] * lcr[0].Stri[st1];
         comm_sums(&ct.tcharge, &ione, COMM_EN2);

        if (pct.gridpe == 0)
            printf ("\n total charge, %18.12f %18.12f  Fermi energy %16.10f %16.10f",
                    ct.nel, ct.tcharge, bias1, density);


        if (fabs(density) < 1.0e-10 && fermi_step == 0)
        {
                bias1a = bias1;
                tchargea = ct.tcharge;
                bias1 += copysign(0.05, ct.nel- ct.tcharge);
        }
        else
        {
            if(fermi_step != 0) density = (ct.tcharge-tchargea)/(bias1 - bias1a);
            bias1a = bias1;
            bias1 += (ct.nel - ct.tcharge)/density;
            tchargea = ct.tcharge;
        }

    }

    if (pct.gridpe == 0)
        printf("\nFERMI ENERGY = %15.8f\n", bias1);

    set_energy_weight (lcr[1].ene, lcr[1].weight, bias1, &lcr[1].nenergy);

    sigma_all_energy_point (sigma_all, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2]);

    lcr[1].EF_new = bias1;
    lcr[2].EF_new = bias1;
    lcr[1].EF_old = bias1;
    lcr[2].EF_old = bias1;

}
