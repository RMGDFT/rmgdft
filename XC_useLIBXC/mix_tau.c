#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "common_prototypes.h"
#include "input.h"


void mix_tau (double * new_tau, double * tau, int length, int length_x, int length_y, int length_z)
{
    double t1;
    int step, inc = 1;
    static double **rhohist=NULL, **residhist=NULL;

    /*Linear Mixing*/
    if (verify("charge_mixing_type","Linear"))
    {
	
	/* Scale old charge density first*/
	t1 = 1.0 - ct.mix;
	QMD_dscal (length, t1, tau, inc); 

	/*Add the new density*/
	QMD_daxpy (length, ct.mix, new_tau, inc, tau, inc);
    }
    //else {
//	if (verify("charge_mixing_type","Pulay"))
//	{
//	    step = ct.scf_steps;

//	    if (ct.charge_pulay_refresh)
//		step = ct.scf_steps % ct.charge_pulay_refresh;

	    /*Use pulay mixing, result will be in rho*/
//	    pulay_rho(step, length, length_x, length_y, length_z, new_tau, tau, ct.charge_pulay_order, &rhohist, &residhist, ct.charge_pulay_special_metrics, ct.charge_pulay_special_metrics_weight);
	    
//	}
	    
  //  }



}                               /* end mix_tau */

/******/
