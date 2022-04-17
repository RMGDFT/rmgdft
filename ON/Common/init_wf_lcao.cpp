/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
     	Just generates a random start.
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "transition.h"
#include "init_var.h"




#if RMG_FAST_MATH
   #include "fastonebigheader.h"
#else
   #define fasterlog log 
   #define fasterexp exp
#endif

static void atomic_wave_to_orbital(STATE *st, double *phi, SPECIES *sp, int ip, int l, int m, int l_extra);
static void get_one_orbital(STATE *states, int state, double *phi);

void init_wf_lcao(STATE * states)
{
    double *phi = new double[ct.max_orbit_size];

    if (pct.gridpe == 0)
        printf(" LCAO initial wavefunction \n");
    MPI_Barrier(pct.img_comm);

    if(ct.LocalizedOrbitalLayout == LO_projection)
    {
        for (int st = 0; st < LocalOrbital->num_thispe; st++)
        {
            int st_glob = LocalOrbital->index_proj_to_global[st];
            if(st_glob < 0) continue;
            get_one_orbital(states, st_glob, phi);
            
            LocalOrbital->AssignOrbital(st, phi);
        }
    }
    else
    {
        for (int state = ct.state_begin; state < ct.state_end; state++)
        {
            get_one_orbital(states, state, phi);
            for(int idx = 0; idx < states[state].size; idx++)
                states[state].psiR[idx] = phi[idx];
        }
    }
    delete [] phi;
    if (pct.gridpe == 0)
        printf(" LCAO initial wavefunction  down\n");
}

static void get_one_orbital(STATE *states, int state, double *phi)
{
    int ion, species, ist;
    STATE *st;
    SPECIES *sp;
    int state_count, ip, l, m;


    ion = states[state].atom_index;
    species = Atoms[ion].species;
    sp = &Species[species];
    ist = states[state].atomic_orbital_index;
    st = &states[state];

    //count how many atomic wave functions we have
    // if the number of atomic wave functions is not enough, 
   // use the same radial function with l+1;

    state_count= 0;
    int l_extra = 0;
    while(ist >= state_count)
    {
        for (ip = 0; ip < sp->num_atomic_waves; ip++)
        {
            l = sp->atomic_wave_l[ip]+l_extra;
            for (m=0; m < 2*l+1; m++)
            {
                if(ist == state_count)
                {
                    atomic_wave_to_orbital(st, phi, sp, ip, l, m, l_extra);
                }
                state_count++;
            }
        }
        l_extra++;
    }

}                               


static void atomic_wave_to_orbital(STATE *st, double *phi, SPECIES *sp, int ip, int l, int m, int l_extra)
{

    int idx, ix, iy, iz;
    int ixx, iyy, izz;

    double crds[3], hx, hy, hz, x, y, z;
    int  i_r;
    double a,b,c, r, vector[3];

    double r1, r2, fradius, coef1, coef2;
    double dr = sp->r[1] - sp->r[0];

    hx = get_hxgrid() * get_xside();
    hy = get_hygrid() * get_yside();
    hz = get_hzgrid() * get_zside();

    crds[0] = st->crds[0];
    crds[1] = st->crds[1];
    crds[2] = st->crds[2];
    ixx = st->ixmax - st->ixmin + 1;
    iyy = st->iymax - st->iymin + 1;
    izz = st->izmax - st->izmin + 1;



    b = log((sp->r[2] - sp->r[1])/(sp->r[1] - sp->r[0]));
    c = (sp->r[0] * exp(b) - sp->r[1])/(1.0 -exp(b) );
    a = sp->r[0] + c;


    for(ix = 0; ix < ixx; ix++ )
    {
        x = (ix + st->ixmin) *hx - crds[0];
        for(iy = 0; iy < iyy; iy++ )
        {
            y = (iy + st->iymin) *hy - crds[1];
            for(iz = 0; iz < izz; iz++ )
            {
                z = (iz + st->izmin) *hz - crds[2];
                r = sqrt(x*x + y*y + z*z);
                idx = ix * iyy * izz + iy * izz + iz;

                vector[0] = x;
                vector[1] = y;
                vector[2] = z;

                if(r > sp->aradius[ip])  
                {
                    phi[idx] = 0.0;
                    continue;
                }

                if(r < sp->r[0])
                {
                    fradius = sp->atomic_wave[ip][0];
                    vector[0] = sp->r[0];
                    vector[1] = sp->r[0];
                    vector[2] = sp->r[0];
                }
                else
                {
                    if(sp->gtype)
                    {
                        i_r = (int)(r/dr);
                        if(i_r >= sp->rg_points) i_r = sp->rg_points - 1;
                    }
                    else
                    {
                        i_r = (int)(fasterlog ( (r+c)/a) /b);
                    }


                    r1 = sp->r[i_r];
                    r2 = sp->r[i_r+1];
                    coef1 = (r2-r)/(r2-r1);
                    coef2 = (r-r1)/(r2-r1);

                    fradius = coef1 * sp->atomic_wave[ip][i_r]
                        + coef2 * sp->atomic_wave[ip][i_r+1];
                }

                phi[idx] = fradius * Ylm(l, m, vector) *std::pow(r, l_extra);
                if(l_extra > 0) phi[idx] *= exp(-r*6.0);


            }
        }
    }

}
