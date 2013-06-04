/************************** SVN Revision Information **************************
 **    $Id: get_wave 2013-06-02 15:38:05Z BTAN $    **
******************************************************************************/
/*

Get a particular wave st and store it in wave_global

*/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "md.h"

#include "my_scalapack.h"


void get_wave(int st, STATE * states, double *wave_global)
{
	int  j, idx, ione = 1;
	REAL t2;
	REAL *wave_temp;
	register double tcharge;
        double charge_from_wave;
	int mxllda;
	mxllda = MXLLDA;
	double one = 1.0, zero = 0.0;

	/* for parallel libraries */
	int n2 = (ct.state_end-ct.state_begin)* ct.num_states;
	REAL *psi1,  scale;
	int i, st1, state_per_proc;
	state_per_proc = ct.state_per_proc + 2;

	int IA=1,JA=1,IB=1,JB=1, numst = ct.num_states;
	int st11;

	my_malloc_init( wave_temp, P0_BASIS, REAL );
	if (pct.gridpe == 0)
		printf(" Compute %d th wave\n", st);

	printf("print out zz_dis right before transpose to cc :\n");

/*
        for (j=0;  j< MXLLDA * MXLCOL; j++)
	{
		printf("proc %d:   zz_dis[%d] = %f\n", pct.gridpe, j, zz_dis[j]);
	}
*/        
        PDTRAN(&numst, &numst, &one, zz_dis, &IA, &JA, pct.desca, &zero, cc_dis, &IB, &JB, pct.desca);//transpose zz_dis to cc_dis

	printf("print out cc_dis right after transpose to cc :\n");

/*
        for (j=0;  j< MXLLDA * MXLCOL; j++)
	{
		printf("proc %d:   cc_dis[%d] = %f\n", pct.gridpe, j, cc_dis[j]);
	}
*/        

	Cpdgemr2d(numst, numst, cc_dis, IA, JA, pct.desca, coefficient_matrix_row, IB, JB,
			pct.descb, pct.desca[1]);


/*
	printf("print coefficient_matrix_row right after Cpdgemr2d in proc %d:\n", pct.gridpe);
        for (j=0;  j< numst*(ct.state_end - ct.state_begin); j++)
	{
		printf("proc %d:  coefficient[%d] = %f\n",pct.gridpe, j, coefficient_matrix_row[j]);
	}

*/        

	for (idx = 0; idx < NX_GRID * NY_GRID * NZ_GRID; idx++)
		wave_global[idx] = 0.;

	for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
	{ 
                int i;
		st11 = st1 - ct.state_begin;
		scale =  coefficient_matrix_row[st11 * ct.num_states + st]; //column st store the coefficients for st wave
		psi1 = states[st1].psiR;
/*
 *  if you dont trust add_orbit_to_wave, 
 *  you can test it with the following trick to get the summation of orbitals as well
 *
		REAL * psi2;
		my_malloc(psi2, states[st1].size, REAL);
		for(i=0; i< states[st1].size; i++)
			psi2[i] = 1.0;
		density_orbit_X_orbit(st1, st1, scale, psi1, psi2, wave_global, 0, states);
                my_free(psi2);         
*/       
 
		add_orbit_to_wave(st1, scale, psi1, wave_global, states);

	}

	my_barrier();

	idx = NX_GRID * NY_GRID * NZ_GRID;

	global_sums(wave_global, &idx, pct.grid_comm);

        global_to_distribute(wave_global, wave_temp);

/*
 * for debug purpose, let print out the total charge from this wave
 */
	tcharge = 0.0;
	charge_from_wave = 0.0;

	for (idx = 0; idx < P0_BASIS; idx++)
		tcharge += wave_temp[idx] * wave_temp[idx];

	charge_from_wave = real_sum_all(tcharge, pct.grid_comm);


	charge_from_wave *= ct.vel;

	if (pct.gridpe == 0)
		printf("\n total charge from %d wave = %f  with ct.vel = %f \n", st, charge_from_wave, ct.vel);

        my_free(wave_temp);
}
