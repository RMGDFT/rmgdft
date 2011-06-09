#include "main.h"
#include <stdio.h>
#include <float.h>
#include <math.h>


void get_psi_overlaps(REAL *psi_array, REAL *overlap, int numst, int maxst, int numpt, int maxpt)
{

   REAL rzero = 0.0;
   REAL rone = 1.0;
   char *transt = "t";
   char *uplo = "l";

   ssyrk( uplo, transt, &numst, &numpt, &rone, psi_array, &numpt, &rzero, overlap, &numst);

}


