#include "Bdef.h"

void BI_ivmcopy(Int m, Int n, Int *A, Int lda, Int *buff)
/*
 *  performs an vector to matrix copy (unpack) for the data type Int
 */
{
   Int i, j;

   if ( (m == lda) || (n == 1) )
   {
      m = n * m;
      for (i=0; i < m; i++) A[i] = buff[i];
   }
   else if (m == 1)
   {
      for (j=0; j < n; j++) A[j*lda] = buff[j];
   }
   else
   {
      for (j=0; j< n; j++)
      {
         for (i=0; i < m; i++) A[i] = buff[i];
         A += lda;
         buff += m;
      }
   }
}
