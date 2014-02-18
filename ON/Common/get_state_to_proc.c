/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


get_state_to_proc.c


*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main_on.h"


void get_state_to_proc(STATE * states)
{
    int st;

    /* determing how many orbitals one processor should have */
    ct.state_per_proc = (ct.num_states + NPES -1) /NPES;


    my_malloc(state_begin, NPES, int);
    my_malloc(state_end, NPES, int);
    for (st = 0; st < ct.num_states; st++)
    {
        states[st].pe = st/ct.state_per_proc;
        states[st].index = st;
    }

    ct.state_begin = ct.state_per_proc * pct.gridpe;
    ct.state_end = ct.state_begin + ct.state_per_proc;

    if(ct.state_begin > ct.num_states ) ct.state_begin = ct.num_states;
    if(ct.state_end > ct.num_states ) ct.state_end = ct.num_states;

    for(st = 0; st < NPES; st++)
    {

        state_begin[st] = ct.state_per_proc * st;
        state_end[st] = state_begin[st] + ct.state_per_proc ;

        if(state_begin[st] > ct.num_states ) state_begin[st] = ct.num_states;
        if(state_end[st] > ct.num_states ) state_end[st] = ct.num_states;
    }

    st = (ct.state_end - ct.state_begin) * ct.num_states;
    if(st == 0 ) st = ct.num_states;
    my_malloc(Hij_00, st, double);
    my_malloc(Bij_00, st, double);
    my_malloc(work_matrix_row, st, double);
    my_malloc(coefficient_matrix_row, st, double);


}
