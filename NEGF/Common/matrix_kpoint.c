/*
 **    $Id$    **
 ******************************************************************************/


/*
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>




void matrix_kpoint (int size, complex double *matrix_k, double *matrix_yz, double kvecy, double kvecz)
{
/*===========   all matrix elements are splitted according to ========
 * (1) if two overlaping orbitals are in the same cell, the value stored
 *      in Matrix_yz[0+ij],  ij is the index of the matrix element.
 * (2) if two orbitals overlaps after one translate by unit vector in y
 *     and/or z, the values are stored in Matrix_yz[K+ij]
 *     K = size of matrix  * m 
 *     m = 1 for +y translation
 *     m = 2 for -y translation
 *     m = 3 for +z translation
 *     m = 4 for -z translation
 *     m = 5 for +y+z translation
 *     m = 6 for +y-z translation
 *     m = 7 for -y+z translation
 *     m = 8 for -y-z translation
 *     in return matrix_k[ij] = matrix_yz[0+ij] +exp(i*kvecy)
 *     matrix_yz[K +ij] + exp(-i*kvecy) * matrix_yz[2K + ij] ... for all
 *     translatrions
 */
    

    int idx, m;
    double complex  ctem[9];

    ctem[0] = 1.0;
    ctem[1] = cexp(+I * kvecy);
    ctem[2] = cexp(-I * kvecy);
    ctem[3] = cexp(+I * kvecz);
    ctem[4] = cexp(-I * kvecz);
    ctem[5] = cexp(+I * kvecy + I * kvecz);
    ctem[6] = cexp(+I * kvecy - I * kvecz);
    ctem[7] = cexp(-I * kvecy + I * kvecz);
    ctem[8] = cexp(-I * kvecy - I * kvecz);

    for(idx = 0; idx <size; idx++)
        matrix_k[idx] = 0.0;
    for (m = 0; m < 9; m++)
    {
        for(idx = 0; idx <size; idx++)
            matrix_k[idx] += ctem[m] * matrix_yz[m * size +idx];
    }

}
