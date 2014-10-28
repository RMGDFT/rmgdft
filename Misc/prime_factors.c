/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/


#include <math.h>
#include <stdio.h>
#include "main.h"

// find a list of prime factors of a positive number n
int prime_factors(int n, int *factors)
{

    int n_factors, i, j;

    if(n <= 0) 
    {
        n_factors = 1;
        factors[0] = n;
        printf("\n WARNING: trying to find a factorization of a negative number");  
    }
    n_factors = 0;

//  find how many of power 2
    while (n%2 == 0)
    {
        factors[n_factors] = 2;
        n_factors++;
        n = n/2;
    }

//  find how many powe of i
    for(i=3; i <= sqrt(n); i=i+2)
    {
        while (n%i == 0)
        {
            factors[n_factors] = i;
            n_factors++;
            n = n/i;
        }
    }

    // if n is a prime number, so it is
    if(n_factors == 0)
    {
        n_factors = 1;
        factors[0] = n;
    }

    return n_factors;
}

