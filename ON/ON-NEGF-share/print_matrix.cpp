/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include "main.h"
#include "prototypes_on.h" 

#define MAX_COL		10


/************************************************************/
void print_matrix(double *b, int n, int ldb)
{
    int i, j, rows, cols, idx;

    rows = n / MAX_COL;
    cols = n % MAX_COL;

    printf("MATRIX \n");
    for (i = 0; i < n; i++)
    {
        for (idx = 0; idx < rows; idx++)
        {
            for (j = 0; j < MAX_COL; j++)
            {
                printf("%15.8e \t", b[i * ldb + j + idx * MAX_COL]);
            }
            if (idx + 1 < rows || cols > 0)
                printf("\n");
        }
        for (j = 0; j < cols; j++)
        {
            printf("%15.8e \t", b[i * ldb + j + rows * MAX_COL]);
        }
        printf("\n");
    }
    printf("\n");

    fflush(NULL);

}

/************************************************************/

void print_matrix_matlab(char *name, int n, double *a)
{
    int i, j;



    printf("\n %s = ...\n", name);
    printf(" [");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf(" %e ", a[n * i + j]);
            if (j < n - 1)
            {
                printf(",");
                if (!((j + 1) % 7))
                    printf("... \n");
            }
        }
        if (i < n - 1)
            printf(" ;... \n");
    }
    printf(" ]\n");

    fflush(NULL);

}
