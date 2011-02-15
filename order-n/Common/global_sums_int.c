/************************** SVN Revision Information **************************
 **    $Id: global_sums_int.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*

                global_sums.c

  Sums an array over all processors. 
  For serial machines it just returns.



*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


#if MPI


void global_sums(REAL * vect, int *length)
{
    int sizr, steps, blocks, newsize;
    REAL *rptr, *rptr1;
    REAL rptr2[100];
    int ione = 1;

    /* Check for small vector case and handle on stack */
    if (*length < 100)
    {

        sizr = *length;
        dcopy(&sizr, vect, &ione, rptr2, &ione);
        MPI_Allreduce(rptr2, vect, sizr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return;

    }


    my_malloc_init( rptr, MAX_PWRK, REAL );
    newsize = MAX_PWRK;
    blocks = *length / newsize;
    sizr = (*length % newsize);

    rptr1 = vect;

    for (steps = 0; steps < blocks; steps++)
    {

        dcopy(&newsize, rptr1, &ione, rptr, &ione);
        MPI_Allreduce(rptr, rptr1, newsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        rptr1 += newsize;

    }                           /* end for */


    if (sizr)
    {

        dcopy(&sizr, rptr1, &ione, rptr, &ione);
        MPI_Allreduce(rptr, rptr1, sizr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    }                           /* end if */


    my_free(rptr);

}                               /* end global_sums */


void global_sums_int(int *vect, int *length)
{

    int sizr, steps, blocks, newsize;
    int idx;
    int *rptr, *rptr1;
    int rptr2[100];

    my_barrier();
    /* Check for small vector case and handle on stack */
    if (*length < 100)
    {

        sizr = *length;
        my_malloc_init( rptr, sizr, int );
        for (idx = 0; idx < sizr; idx++)
            rptr[idx] = vect[idx];

        MPI_Allreduce(rptr, vect, sizr, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

        my_free(rptr);
        return;

    }                           /* end if */


    my_malloc_init( rptr, MAX_PWRK, int );
    newsize = MAX_PWRK;
    blocks = *length / newsize;
    sizr = (*length % newsize);

    rptr1 = vect;

    for (steps = 0; steps < blocks; steps++)
    {

        for (idx = 0; idx < newsize; idx++)
            rptr[idx] = rptr1[idx];
        MPI_Allreduce(rptr, rptr1, newsize, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

        rptr1 += newsize;

    }                           /* end for */


    if (sizr)
    {
        for (idx = 0; idx < sizr; idx++)
            rptr[idx] = rptr1[idx];

        MPI_Allreduce(rptr, rptr1, sizr, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

    }                           /* end if */


    my_free(rptr);


}                               /* end global_sums_int */


#endif /*  endif MPI  */


#if SERIAL

void global_sums(double *vect, int *length)
{

    return;

}                               /* end global_sums */

void global_sums_int(int *vect, int *length)
{
    return;
}

#endif
