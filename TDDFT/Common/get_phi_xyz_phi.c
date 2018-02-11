/************************** SVN Revision Information **************************
 **    $Id: 
******************************************************************************/

/* calculating Hamiltonian matrix Hij and 
 * overlap matrix matB togather
 */

 
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "Scalapack.h"



void get_phi_xyz_phi(STATE * states, double *Xij, double *Yij, double *Zij)
{
    int idx, st1, st2, idx1, idx2;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    double tem, vel;
    int ixx, iyy, izz;
    double *Xij_00, *Yij_00, *Zij_00;

    int IA=1, JA=1, IB=1, JB=1, numst = ct.num_states;

    idx1 = ct.num_states * ct.num_states;

    maxst = ct.num_states;


    for (st1 = 0; st1 < (ct.state_end-ct.state_begin) * ct.num_states; st1++)
    {
        Xij[st1] = 0.;
        Yij[st1] = 0.;
        Zij[st1] = 0.;
    }



    /* calculate the < states.psiR |x,y,z| states1.psiR> for electric
 * field kick in TDDFT  */


    //orbit_xyz_orbit(states, Xij_00, Yij_00, Zij_00);
    orbit_xyz_orbit(states, Xij, Yij, Zij);


    //my_barrier();

    n2 = (ct.state_end-ct.state_begin) * ct.num_states;
    vel = get_vel();
    dscal (&n2, &vel, Xij, &ione);
    dscal (&n2, &vel, Yij, &ione);
    dscal (&n2, &vel, Zij, &ione);

    if (pct.gridpe == 0)
    {
        printf(" matrix Xij \n");
        print_matrix(Xij, 5, maxst);
        print_matrix(Yij, 5, maxst);
        print_matrix(Zij, 5, maxst);
    }

/*
    Cpdgemr2d(numst, numst, Xij_00, IA, JA, pct.descb, Xij, IB, JB,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Yij_00, IA, JA, pct.descb, Yij, IB, JB,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Zij_00, IA, JA, pct.descb, Zij, IB, JB,
            pct.desca, pct.desca[1]);
*/


}
