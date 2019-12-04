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
#include "init_var.h"




#if RMG_FAST_MATH
   #include "fastonebigheader.h"
#else
   #define fasterlog log 
   #define fasterexp exp
#endif

static void atomic_wave_to_orbital(STATE *st, double *phi, SPECIES *sp, int ip, int l, int m);
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
    int idx;
    int ion, species, ist;
    STATE *st;
    SPECIES *sp;
    int state_count, ip, l, m, ix, iy, iz, ixx, iyy, izz;
    int idx1, idx2, idx3, idx4, idx5, idx6;
    long idum;


    ion = states[state].atom_index;
    species = Atoms[ion].species;
    sp = &Species[species];
    ist = states[state].atomic_orbital_index;
    st = &states[state];

    //count how many atomic wave functions we have

    state_count= 0;
    for (ip = 0; ip < sp->num_atomic_waves; ip++)
    {
        l = sp->atomic_wave_l[ip];
        for (m=0; m < 2*l+1; m++)
        {
            if(ist == state_count)
            {
                atomic_wave_to_orbital(st, phi, sp, ip, l, m);
            }
            state_count++;
        }
    }


    // random start for this orbitals
    if(ist >= state_count) 
    {

        idum = 3356 + state;  /*  seeds are different for different orbitals  */
        rand0(&idum);
        ixx = states[state].orbit_nx;
        iyy = states[state].orbit_ny;
        izz = states[state].orbit_nz;
        for(ix = 0; ix < ixx; ix++)
            for(iy = 0; iy < iyy; iy++)
                for(iz = 0; iz < izz; iz++)
                {
                    idx = ix * iyy * izz + iy * izz + iz;
                    phi[idx] = 0.0;
                }

        for(ix = ixx/2 -3; ix < ixx/2+3; ix++)
            for(iy = iyy/2 -3; iy < iyy/2+3; iy++)
                for(iz = izz/2 -3; iz < izz/2+3; iz++)
                {
                    idx = ix * iyy * izz + iy * izz + iz;
                    phi[idx] = rand0(&idum);
                }

        for(ix = ixx/2 -4; ix < ixx/2+4; ix++)
            for(iy = iyy/2 -4; iy < iyy/2+4; iy++)
                for(iz = izz/2 -4; iz < izz/2+4; iz++)
                {
                    idx = ix * iyy * izz + iy * izz + iz;
                    idx1 = (ix-1) *iyy * izz + (iy+0) * izz + iz +0;
                    idx2 = (ix+1) *iyy * izz + (iy+0) * izz + iz +0;
                    idx3 = (ix+0) *iyy * izz + (iy-1) * izz + iz +0;
                    idx4 = (ix+0) *iyy * izz + (iy+0) * izz + iz +0;
                    idx5 = (ix+0) *iyy * izz + (iy+0) * izz + iz -1;
                    idx6 = (ix+0) *iyy * izz + (iy+0) * izz + iz +1;

                    phi[idx] += (phi[idx1] +phi[idx2] +phi[idx3]
                            +phi[idx4] +phi[idx5] +phi[idx6])/6.0 ;
                }


    }


}                               


static void atomic_wave_to_orbital(STATE *st, double *phi, SPECIES *sp, int ip, int l, int m)
{

    int idx, ix, iy, iz;
    int ixx, iyy, izz;

    double crds[3], hx, hy, hz, x, y, z;
    int  i_r, yindex;
    double a,b,c, r, vector[3];

    double r1, r2, fradius, coef1, coef2;
    double dr = sp->r[1] - sp->r[0];



    yindex = l*l + m;
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
                    fradius = sp->atomic_wave[l][0];
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

                phi[idx] = fradius * ylm(yindex, vector);


            }
        }
    }

}
