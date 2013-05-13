#include "main.h"
#include <stdio.h>
#include <float.h>
#include <math.h>


void get_psi_overlaps(rmg_double_t *psi_array, rmg_double_t *overlap, int numst, int maxst, int numpt, int maxpt)
{

   rmg_double_t rzero = 0.0;
   rmg_double_t rone = 1.0;
   char *transt = "t";
   char *uplo = "l";

   ssyrk( uplo, transt, &numst, &numpt, &rone, psi_array, &numpt, &rzero, overlap, &numst);

}


