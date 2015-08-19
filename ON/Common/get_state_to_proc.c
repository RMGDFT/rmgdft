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
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"

void get_state_to_proc(STATE * states)
{
    int st, pe, proc_extra;

    assert(ct.num_states >= pct.grid_npes);
    /* determing how many orbitals one processor should have */
    ct.state_per_proc = (ct.num_states + pct.grid_npes -1) /pct.grid_npes;


    my_malloc(state_begin, pct.grid_npes, int);
    my_malloc(state_end, pct.grid_npes, int);
    for (st = 0; st < ct.num_states; st++)
    {
        states[st].index = st;
    }
#if ELEMENTAL_LIBS
    //  first proc_extra processor will have one more orbitals than other processors.
    proc_extra = ct.num_states%pct.grid_npes;
    if (proc_extra == 0)
    {
        for(pe = 0; pe < pct.grid_npes; pe++)
        {
            state_begin[pe] = ct.state_per_proc * pe;
            state_end[pe] = ct.state_per_proc * (pe+1);
        }

    }
    else
    {
        for(pe = 0; pe < proc_extra; pe++)
        {
            state_begin[pe] = ct.state_per_proc * pe;
            state_end[pe] = ct.state_per_proc * (pe+1);
        }

        for(pe = proc_extra; pe < pct.grid_npes; pe++)
        {
            state_begin[pe] = state_end[pe-1];
            state_end[pe] = state_begin[pe] + ct.num_states/pct.grid_npes;
        }
    }

    ct.state_begin = state_begin[pct.gridpe];
    ct.state_end = state_end[pct.gridpe];
    if(ct.state_end > ct.num_states)
    {
        dprintf("\n something wroth in get_state_to_proc  %d %d", ct.state_end, pct.gridpe);
        exit(0);
    }

#else

    ct.state_begin = ct.state_per_proc * pct.gridpe;
    ct.state_end = ct.state_begin + ct.state_per_proc;

    if(ct.state_begin > ct.num_states ) ct.state_begin = ct.num_states;
    if(ct.state_end > ct.num_states ) ct.state_end = ct.num_states;

    for(st = 0; st < pct.grid_npes; st++)
    {

        state_begin[st] = ct.state_per_proc * st;
        state_end[st] = state_begin[st] + ct.state_per_proc ;

        if(state_begin[st] > ct.num_states ) state_begin[st] = ct.num_states;
        if(state_end[st] > ct.num_states ) state_end[st] = ct.num_states;
    }
#endif

    for(pe = 0; pe < pct.grid_npes; pe++)
        for (st = state_begin[pe]; st < state_end[pe]; st++)
            states[st].pe = pe;

    st = (ct.state_end - ct.state_begin) * ct.num_states;
    if(st == 0 ) st = 1;
    my_malloc(Hij_00, st, double);
    my_malloc(Bij_00, st, double);
    my_malloc(theta, st, double);
    my_malloc(work_matrix_row, st, double);
    my_malloc(coefficient_matrix_row, st, double);


}
