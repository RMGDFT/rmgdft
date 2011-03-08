/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void projection (double *data, int x, int y, int z, char direction, char *tag)
{
    double p[10000];
    int i, j, k;
    double t;

    for (i = 0; i < 10000; i++)
        p[i] = 0;
    t = 0;

    if (direction == 'x' || direction == 'X')
    {
        for (i = 0; i < x; i++)
        {
            for (j = 0; j < y; j++)
            {
                for (k = 0; k < z; k++)
                    p[i] += data[i * y * z + j * z + k];
            }
        }
        for (i = 0; i < x; i++)
        {
            printf ("\n");
            printf (tag);
            printf ("%d:  %f ", i, p[i]);
            t += p[i];
        }
        fflush (NULL);
    }

    if (direction == 'y' || direction == 'Y')
    {
        for (i = 0; i < x; i++)
        {
            for (j = 0; j < y; j++)
            {
                for (k = 0; k < z; k++)
                {
                    p[j] += data[i * y * z + j * z + k];
                }
            }
        }
        for (i = 0; i < y; i++)
        {
            printf ("\n");
            printf (tag);
            printf (" %f ", p[i]);
        }
        fflush (NULL);
    }

    if (direction == 'z' || direction == 'Z')
    {
        for (i = 0; i < x; i++)
        {
            for (j = 0; j < y; j++)
            {
                for (k = 0; k < z; k++)
                {
                    p[k] += data[i * y * z + j * z + k];
                }
            }
        }
        for (i = 0; i < z; i++)
        {
            printf ("\n");
            printf (tag);
            printf (" %f ", p[i]);
        }
        fflush (NULL);
    }

}
