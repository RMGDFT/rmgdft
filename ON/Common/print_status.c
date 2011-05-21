/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void print_status(STATE * states, REAL * vh, REAL * vxc, REAL * vnuc, REAL * vh_old, char *msg)
{
    int i, j;
    REAL total;

    printf("\n print_status ----%s ", msg);

    total = 0;
    for (i = 0; i < P0_BASIS; i++)
        total += vh[i];
    printf("PE: %d, vh total: %22.16f \n", pct.gridpe, total);

    total = 0;
    for (i = 0; i < P0_BASIS; i++)
        total += vxc[i];
    printf("PE: %d, vxc total: %22.16f \n", pct.gridpe, total);

    total = 0;
    for (i = 0; i < P0_BASIS; i++)
        total += vnuc[i];
    printf("PE: %d, vnuc total: %22.16f \n", pct.gridpe, total);

    total = 0;
    for (i = 0; i < P0_BASIS; i++)
        total += vh_old[i];
    printf("PE: %d, vh_old total: %22.16f \n", pct.gridpe, total);

    for (i = ct.state_begin; i < ct.state_end; i++)
    {

        total = 0;
        for (j = 0; j < states[i].size; j++)
        {
            total += states[i].psiR[j];
        }
        printf("PE: %d, state: %d, orbital total: %22.16f \n", pct.gridpe, i, total);

    }
    fflush(NULL);

}

void print_state_projections(STATE * states, char direction)
{

    char newname[100];
    int st;


    for (st = ct.state_begin; st < ct.state_end; st++)
    {
        sprintf(newname, "PE%d.STATE%d:", pct.gridpe, st);
/*        projection(states[st].psiR, states[st].orbit_nx, states[st].orbit_ny,
 *                  states[st].orbit_nz, direction, newname);
*/
    }
    fflush(NULL);

}

void print_state_sum(STATE * states)
{
    int st;
    int i;
    REAL temp;

    for (st = ct.state_begin; st < ct.state_end; st++)
    {
        temp = 0.0;
        for (i = 0; i < states[st].size; i++)
            temp += states[st].psiR[i];
        printf("\n PE: %d STATE: %d ---sum---: %22.16f", pct.gridpe, st, temp);
    }
    fflush(NULL);
}

void print_state(STATE * state)
{
    int st;
    int i;
    REAL temp;

    temp = 0.0;
    for (i = 0; i < state->size; i++)
        temp += state->psiR[i];
    printf("\n PE: %d.STATE: %d---sum---: %22.16f", pct.gridpe, temp);
    fflush(NULL);
}

void print_sum(int size, double *data, char *msg)
{
    double sum = 0.0;
    int idx;

    for (idx = 0; idx < size; idx++)
        sum += data[idx];
    printf("\n %s  ----  %22.14f  \n", msg, sum);
    fflush(NULL);

}

void print_data(int size, double *data)
{
    int idx;

    for (idx = 0; idx < size; idx++)
    {
        printf("\n %d    %22.14f  ", idx, data[idx]);
    }
}


void print_sum_square(int size, double *data, char *msg)
{
    double sum = 0.0;
    int idx;

    for (idx = 0; idx < size; idx++)
        sum += data[idx] * data[idx];
    printf("\n %s  ----  %22.16f", msg, sum);

}

void print_states_dot_product(STATE * states)
{
    int st;
    int i;
    REAL temp;
    for (st = ct.state_begin; st < ct.state_end; st++)
    {
        temp = 0.0;
        for (i = 0; i < states[st].size; i++)
            temp += states[st].psiR[i] * states[st].psiR[i];
        printf("\n PE: %d.STATE: %d dot_product : %22.16f", pct.gridpe, st, temp);
    }
    fflush(NULL);
}

void print_sum_idx(int size, double *data, char *msg)
{
    double sum = 0.0;
    int idx;

    for (idx = 0; idx < size; idx++)
        sum += data[idx] * (idx - size / 2);
    printf("\n %s  ----  %22.16f", msg, sum);

}

void print_orbit_centers(STATE * states)
{
    REAL x, y, z;
    int st;

    for (st = ct.state_begin; st < ct.state_end; st++)
    {
        get_orbit_center(&states[st], &x, &y, &z);
        printf
            (" STATE: %d--center: %f, %f, %f; centroid: %f, %f, %f; difference: %f, %f, %f \n",
             st, states[st].crds[0], states[st].crds[1], states[st].crds[2],
             x, y, z, x - states[st].crds[0], y - states[st].crds[1], z - states[st].crds[2]);
    }
}
