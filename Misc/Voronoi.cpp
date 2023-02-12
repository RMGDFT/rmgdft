#include <iostream>     
#include <algorithm>    
#include <cfloat>
#include <math.h>       
#include <mpi.h>       
#include "RmgException.h"
#include "Voronoi.h"
#include "transition.h"
#include "GlobalSums.h"
#include "Atomic.h"

Voronoi::Voronoi()
{
    int density = Rmg_G->default_FG_RATIO;
    pbasis = Rmg_G->get_P0_BASIS(density);
    grid_to_atom = new int[pbasis];

    int PX0_GRID = Rmg_G->get_PX0_GRID(density);
    int PY0_GRID = Rmg_G->get_PY0_GRID(density);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(density);
    int PX_OFFSET = Rmg_G->get_PX_OFFSET(density);
    int PY_OFFSET = Rmg_G->get_PY_OFFSET(density);
    int PZ_OFFSET = Rmg_G->get_PZ_OFFSET(density);
    double hx = Rmg_G->get_hxgrid(density);
    double hy = Rmg_G->get_hygrid(density);
    double hz = Rmg_G->get_hzgrid(density);
    double xc = PX_OFFSET * hx; 
    double yc = PY_OFFSET * hy;
    double zc = PZ_OFFSET * hz;

    double xtal[3], dist, min_dist, xcry[3];
    for(int ix = 0; ix < PX0_GRID; ix++)
        for(int iy = 0; iy < PY0_GRID; iy++)
            for(int iz = 0; iz < PZ0_GRID; iz++)
            {
                int idx = ix * PY0_GRID * PZ0_GRID + iy * PZ0_GRID + iz;
                xtal[0] = xc + ix * hx;
                xtal[1] = yc + iy * hy;
                xtal[2] = zc + iz * hz;

                min_dist = DBL_MAX;
                for(size_t ion= 0; ion < Atoms.size(); ion++)
                {
                    xcry[0] = xtal[0] - Atoms[ion].xtal[0];
                    xcry[1] = xtal[1] - Atoms[ion].xtal[1];
                    xcry[2] = xtal[2] - Atoms[ion].xtal[2];

                    if(xcry[0] > 0.5) xcry[0] -= 1.;
                    if(xcry[1] > 0.5) xcry[1] -= 1.0;
                    if(xcry[2] > 0.5) xcry[2] -= 1.0;
                    if(xcry[0] < -0.5) xcry[0] += 1.0;
                    if(xcry[1] < -0.5) xcry[1] += 1.0;
                    if(xcry[2] < -0.5) xcry[2] += 1.0;
                    
                    dist = Rmg_L.metric(xcry);
                    if(dist < min_dist)
                    {
                        min_dist = dist;
                        grid_to_atom[idx] = ion;
                    }

                }


            }

    localrho_atomic = new double[Atoms.size()];
    double *atomic_rho = new double[ct.nspin*pbasis];
    InitLocalObject (atomic_rho, pct.localatomicrho, ATOMIC_RHO, false);
    if(ct.AFM) 
    {
        Rmg_Symm->symmetrize_rho_AFM(atomic_rho, &atomic_rho[pbasis]);
    }
    LocalCharge(atomic_rho, localrho_atomic);
    delete [] atomic_rho;

}

Voronoi::~Voronoi()
{
    delete [] grid_to_atom;
    delete [] localrho_atomic;
}

void Voronoi::LocalCharge(double *rho, double *localrho)
{ 
    int density = Rmg_G->default_FG_RATIO;
    int NX_GRID = Rmg_G->get_NX_GRID(density);
    int NY_GRID = Rmg_G->get_NZ_GRID(density);
    int NZ_GRID = Rmg_G->get_NY_GRID(density);
    double vol =  Rmg_L.get_omega()/(NX_GRID * NY_GRID * NZ_GRID);


    for(size_t ion = 0; ion < Atoms.size(); ion++) localrho[ion] = 0.0;
    for(int idx = 0; idx < pbasis; idx++) localrho[grid_to_atom[idx]] += vol * rho[idx];
    if(ct.AFM) 
    {
        for(int idx = 0; idx < pbasis; idx++) localrho[grid_to_atom[idx]] += vol * rho[idx + pbasis];
    }
    GlobalSums(localrho, (int)Atoms.size(), pct.grid_comm);  
    GlobalSums(localrho, (int)Atoms.size(), pct.spin_comm);  
}
