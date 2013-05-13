/* Applies non-local and S to all orbitals at once and stores results in pct.nv and pct.ns */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"



void app_nls_batch (STATE *sp, rmg_double_t *nv, rmg_double_t *ns, rmg_double_t *Bns, rmg_double_t *sintR)
{

   int istate;
   rmg_double_t *tmp_psi;


   tmp_psi = sp->psiR;
   app_nls_allstates (tmp_psi, NULL, nv, NULL, ns, NULL, Bns, NULL, sintR, NULL, sp->kidx);

}

