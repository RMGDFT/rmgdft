/************************** SVN Revision Information **************************
 **    $Id: allocate_masks.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void allocate_masks(STATE * states)
{

    int radius_order[MAX_STATES];
    REAL radius[MAX_STATES];
    int num_masks;
    REAL mask_radius[MAX_STATES];
    int mask_size[MAX_STATES];
    int state;
    int temp_i;
    REAL temp_r;
    char *mask_ptr;
    int i;
    int idx;
    int order;


/*
	for(i=0; i<MAX_STATES; i++)
	{
		radius_order[i] = 0; 
		radius[i] = 0; 
		mask_radius[i] = 0; 
		mask_size[i] = 0; 
	}

	for(state=0; state<ct.num_states; state++)
	{
		radius_order[state] = state; 
		radius[state] = states[state].radius; 
	}

	for(order=1; order<ct.num_states; order++)
	{
		for(state=0; state<ct.num_states-order; state++)
		{
			if (radius[state] < radius[state+1]) 
			{
				temp_i = radius_order[state]; 
				radius_order[state] = radius_order[state+1]; 
				radius_order[state+1] = temp_i; 
				temp_r = radius[state]; 
				radius[state] = radius[state+1]; 
				radius[state+1] = temp_r; 
			}
		}
	}

	num_masks = 1; 
	mask_radius[0] = radius[0]; 
	mask_size[0] = 2 * states[radius_order[0]].size; 
	for(state=0; state<ct.num_states-1; state++) 
		if(radius[state] != radius[state+1])
		{
			mask_radius[num_masks] = radius[state+1]; 
			mask_size[num_masks] = 2 * states[radius_order[state+1]].size; 
			num_masks++; 
		}

	for(i=0; i<num_masks; i++) 
	{
		my_malloc( mask_ptr[i], mask_size[i], char ); 
		for(idx=0; idx<mask_size[i]; idx++) 
			mask_ptr[i][idx] = 0; 
	}

*/

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        my_malloc( mask_ptr, 2 * states[state].size, char );
        states[state].lmask[0] = mask_ptr;
    }
}
