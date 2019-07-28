#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"


void Sgreen_c (double * Htri, double * Stri, std::complex<double> * sigma1, std::complex<double> * sigma2,
                    std::complex<double> ene, std::complex<double> * Green_C, int nC)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   Green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    std::complex<double> *H_tri, *ch00;

    int info;
    int i, j;
    int nmax, *ipiv;
    int ni[MAX_BLOCKS], ntot, ndim;
    int N;
    double tem;


    nmax = nC;


    N = ct.num_blocks;
    for (i = 0; i < N; i++)
    {
        ni[i] = ct.block_dim[i];
    }

    ntot = 0;
    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ntot += ni[i] * ni[i];
        ndim += ni[i];
    }

    if (ndim != nC)
    {
        printf ("\n %d %d ndim  not equaol to nC in Sgreen_c.c", ndim, nC);
        exit (0);
    }


    for (i = 1; i < N; i++)
    {
        ntot += ni[i - 1] * ni[i];
    }


    /* allocate matrix and initialization  */
    my_malloc_init( H_tri, ntot, std::complex<double> );
    my_malloc_init( ch00, nmax * nmax, std::complex<double> );
    my_malloc_init( ipiv, nmax, int );

    for (i = 0; i < ni[0] * ni[0]; i++)
    {
        H_tri[i] = ene * Stri[i] - Htri[i] * Ha_eV - sigma1[i];
    }

    for (i = ni[0] * ni[0]; i < ntot - ni[N - 1] * ni[N - 1]; i++)
    {
        H_tri[i] = ene * Stri[i] - Htri[i] * Ha_eV;
    }

    j = 0;
    for (i = ntot - ni[N - 1] * ni[N - 1]; i < ntot; i++)
    {
        H_tri[i] = ene * Stri[i] - Htri[i] * Ha_eV - sigma2[j];
        j++;
    }

    tri_to_whole_complex( H_tri, ch00, ct.num_blocks, ct.block_dim);

    get_inverse_block (ch00, Green_C, ipiv, nmax);


    my_free( H_tri );
    my_free( ch00 );
    my_free( ipiv );

}
