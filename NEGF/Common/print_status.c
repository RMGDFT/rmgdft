/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"


void print_status (STATE * states, double * vh, double * vxc, double * vnuc, double * vh_old, char *msg)
{
    int i, j;
    double total;

    printf ("\n print_status ----%s ", msg);

    total = 0;
    for (i = 0; i < get_P0_BASIS(); i++)
        total += vh[i];
    printf ("PE: %d, vh total: %22.16f \n", pct.gridpe, total);

    total = 0;
    for (i = 0; i < get_P0_BASIS(); i++)
        total += vxc[i];
    printf ("PE: %d, vxc total: %22.16f \n", pct.gridpe, total);

    total = 0;
    for (i = 0; i < get_P0_BASIS(); i++)
        total += vnuc[i];
    printf ("PE: %d, vnuc total: %22.16f \n", pct.gridpe, total);

    total = 0;
    for (i = 0; i < get_P0_BASIS(); i++)
        total += vh_old[i];
    printf ("PE: %d, vh_old total: %22.16f \n", pct.gridpe, total);

    for (i = ct.state_begin; i < ct.state_end; i++)
    {

        total = 0;
        for (j = 0; j < states[i].size; j++)
        {
            total += states[i].psiR[j];
        }
        printf ("PE: %d, state: %d, orbital total: %22.16f \n", pct.gridpe, i, total);

    }
    fflush (NULL);

}


void print_state_sum (STATE * states)
{
    int st;
    int i;
    double temp;

    for (st = ct.state_begin; st < ct.state_end; st++)
    {
        temp = 0.0;
        for (i = 0; i < states[st].size; i++)
            temp += states[st].psiR[i];
        printf ("\n PE: %d STATE: %d ---sum---: %22.16f", pct.gridpe, st, temp);
    }
    fflush (NULL);
}

void print_state (STATE * state)
{
    int st;
    int i;
    double temp;

    temp = 0.0;
    for (i = 0; i < state->size; i++)
        temp += state->psiR[i];
    printf ("\n PE: %d.STATE: %d---sum---: %22.16f", pct.gridpe, temp);
    fflush (NULL);
}

void print_sum (int size, double *data, char *msg)
{
    double sum = 0.0;
    int idx;

    for (idx = 0; idx < size; idx++)
        sum += data[idx];
    printf ("\n %s  ----  %22.14f  \n", msg, sum);
    fflush (NULL);

}

void print_data (int size, double *data)
{
    int idx;

    for (idx = 0; idx < size; idx++)
    {
        printf ("\n %d    %22.14f  ", idx, data[idx]);
    }
}


void print_sum_square (int size, double *data, char *msg)
{
    double sum = 0.0;
    int idx;

    for (idx = 0; idx < size; idx++)
        sum += data[idx] * data[idx];
    printf ("\n %s  ----  %22.16f", msg, sum);

}

void print_states_dot_product (STATE * states)
{
    int st;
    int i;
    double temp;
    for (st = ct.state_begin; st < ct.state_end; st++)
    {
        temp = 0.0;
        for (i = 0; i < states[st].size; i++)
            temp += states[st].psiR[i] * states[st].psiR[i];
        printf ("\n PE: %d.STATE: %d dot_product : %22.16f", pct.gridpe, st, temp);
    }
    fflush (NULL);
}

void print_sum_idx (int size, double *data, char *msg)
{
    double sum = 0.0;
    int idx;

    for (idx = 0; idx < size; idx++)
        sum += data[idx] * (idx - size / 2);
    printf ("\n %s  ----  %22.16f", msg, sum);

}



void print_data_to_file (int size, double *data, char *filename)
{

    int i;
    FILE *file;

    file = fopen (filename, "a");
    if (file == NULL)
        error_handler ("Unable to open file to write ");

    for (i = 0; i < size; i++)
    {
        fprintf (file, "%22.16f  \n", data[i]);
    }
    fclose (file);

}
