/************************** SVN Revision Information **************************
 **    $Id: init_state_distribute.c 3157 2015-08-18 16:45:37Z luw $    **
 ******************************************************************************/

/*

    read orbitals which has overlap with this processor domain and map it on
*/



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"



void InitStatedistribute ()
{
    int fhand;
    long nbytes;
    char newname[MAX_PATH + 200];
    char msg[200];

    int idx, idx0, idx2, subsystem;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    int incx, incy, ix, iy, iz;

    int size, idx1;
    int NX, NY, NZ;
    int ixmin_old, ixmax_old, iymin_old, iymax_old, izmin_old, izmax_old;
    int st, st1, st2, st_new, max_orbit_nx_ny;
    int item;

    int x_off, y_off, z_off;
    int istart, block_i, st_in_block;

    double *psi_whole; 
    double hx_old, hx_new, hy_old, hy_new;
    double x1_old, x1_new, y1_old, y1_new;


    size = ct.max_orbit_nx * ct.max_orbit_ny * ct.max_orbit_nz;
    psi_whole = new double[size];

    hx_new = get_hxgrid() * get_xside();
    hy_new = get_hygrid() * get_yside();

    /* Wait until everybody gets here */
    my_barrier ();

    /*
       printf ("state_begin, state_end %d %d \n", ct.state_begin, ct.state_end);
       */


    x_off = get_PX_OFFSET();
    y_off = get_PY_OFFSET();
    z_off = get_PZ_OFFSET();

//  local_index = -1 means that this orbital is not on this process
    states_distribute = new STATE[ct.num_states];
    for (st = 0; st < ct.num_states; st++) states_distribute[st].local_index = -1;

    pct.num_local_orbit = 0;

    for (st = 0; st < ct.num_states; st++) 
    {

        ixmin = states[st].ixmax;
        ixmax = states[st].ixmax;

        if(ixmin >= 0 && ixmax < get_NX_GRID())
        {
            if(states[st].ixmax <x_off || states[st].ixmin > x_off + get_PX0_GRID()) continue;
        }

        if(ixmin < 0)
        {
            ixmin += get_NX_GRID();
            if(x_off > states[st].ixmax && x_off + get_PX0_GRID() < ixmin) continue;
        }

        if(ixmax >= get_NX_GRID() ) 
        {
            ixmax -= get_NX_GRID();
            if(x_off > ixmax && x_off + get_PX0_GRID() < ixmin) continue;
        }

        iymin = states[st].iymin;
        iymax = states[st].iymax;

        if(iymin >= 0 && iymax < get_NY_GRID())
        {
            if(states[st].iymax <y_off || states[st].iymin > y_off + get_PY0_GRID()) continue;
        }

        if(iymin < 0)
        {
            iymin += get_NY_GRID();
            if(y_off > states[st].iymax && y_off + get_PY0_GRID() < iymin) continue;
        }

        if(iymax >= get_NY_GRID() ) 
        {
            iymax -= get_NY_GRID();
            if(y_off > iymax && y_off + get_PY0_GRID() < iymin) continue;
        }

        izmin = states[st].izmin;
        izmax = states[st].izmax;
        if(izmin >= 0 && izmax < get_NZ_GRID())
        {
            if(states[st].izmax <z_off || states[st].izmin > z_off + get_PZ0_GRID()) continue;
        }

        if(izmin < 0)
        {
            izmin += get_NZ_GRID();
            if(z_off > states[st].izmax && z_off + get_PZ0_GRID() < izmin) continue;
        }

        if(izmax >= get_NZ_GRID() ) 
        {
            izmax -= get_NZ_GRID();
            if(z_off > izmax && z_off + get_PZ0_GRID() < izmin) continue;
        }


        states_distribute[pct.num_local_orbit].istate = st;
        states_distribute[st].local_index = pct.num_local_orbit;
        states_distribute[pct.num_local_orbit].whichblock = block_i;
        states_distribute[pct.num_local_orbit].istate_in_block = st_in_block;
        pct.num_local_orbit++;
    }


    printf("\n PE %d  num_local_orbit  %d", pct.gridpe, pct.num_local_orbit);
    double *rptr;
    size = pct.num_local_orbit * get_P0_BASIS()+1024;

    rptr = new double[size];

    for (st1 = 0; st1 < pct.num_local_orbit; st1++)
    {
        states_distribute[st1].psiR = rptr;
        rptr += get_P0_BASIS();
    }




    for (st2 = 0; st2 < pct.num_local_orbit; st2++)
    {
        st = states_distribute[st2].istate;

        sprintf (newname, "%s%s%d", ct.infile, ".orbit_", st);
        fhand = open(newname, O_RDWR, S_IREAD | S_IWRITE );
        if (fhand < 0)
        {
            printf ("\n unable to open file  %s, st1 = %d %d", newname, st1, st);
            exit(0);
        }

        idx = states[st].size * (int) sizeof (double);

        /* ====================== Reading orbitals ======================= */                    

        /*
           printf ("reading orbitals for 3rd/4th probe, idx0 =  %d %s \n", idx0, newname);
         */

        nbytes = read (fhand, psi_whole, idx);
        nbytes = read (fhand, &ixmin_old, sizeof (int));
        nbytes = read (fhand, &ixmax_old, sizeof (int));
        nbytes = read (fhand, &iymin_old, sizeof (int));
        nbytes = read (fhand, &iymax_old, sizeof (int));
        if (nbytes != sizeof (int))
        {
            printf ("\n read %d is different from %d for state %d\n", (int) nbytes,
                    (int) sizeof (int), st);
            exit(0);
        }

        close (fhand);

        MapOrbitalToProcess(st2, states, states_distribute, psi_whole);
    }   

}   


