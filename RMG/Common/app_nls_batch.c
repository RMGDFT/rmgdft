/* Applies non-local and S to all orbitals at once and stores results in pct.nv and pct.ns */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"



void app_nls_batch (STATE *sp, REAL *nv, REAL *ns, REAL *sintR)
{

   int istate;
   REAL *tmp_psi;


   tmp_psi = sp->psiR;
   app_nls_allstates (tmp_psi, NULL, nv, NULL, ns, NULL, sintR, NULL, sp->kidx);

}

