/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/


#include <math.h>
#include <stdio.h>

// find a list of prime factors of a positive number n
int prime_factors(int n, int *factors)
{

    int n_factors, i;

    if(n <= 1) 
    {
        n_factors = 1;
        factors[0] = n;
        printf("\n WARNING: trying to find a factorization of a negative number");  
    }
    n_factors = 0;

//  find how many powe of i
    for(i=2; i <= n; i++)
    {
        while (n%i == 0)
        {
            factors[n_factors] = i;
            n_factors++;
            n = n/i;
        }
    }

    return n_factors;
}

