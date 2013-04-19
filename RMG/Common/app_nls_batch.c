/* Applies non-local and S to all orbitals at once and stores results in pct.nv and pct.ns */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"



void app_nls_batch (STATE *sp, REAL *nv, REAL *ns)
{

   int istate;
   REAL *tmp_psi;

 
   for(istate = 0;istate < ct.num_states;istate++) {

       tmp_psi = sp->psiR;
       app_nls (tmp_psi, NULL, &nv[istate * pct.P0_BASIS], NULL, &ns[istate * pct.P0_BASIS], NULL, pct.oldsintR_local, NULL, sp->istate, sp->kidx);
       sp++;

   }

}

