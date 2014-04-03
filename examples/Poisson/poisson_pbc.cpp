/**
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 *   All rights reserved.
 */


#include <math.h>
#include "BaseThread.h"
#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "vhartree.h"
#include "RmgSumAll.h"
#include "RmgTimer.h"
#include "rmg_error.h"
#include "boundary_conditions.h"
#include <mpi.h>

using namespace std;


string header =
"\n"
" 3D parallel poission solver:\n"
"   Usage  mpirun -np 8 ./poisson_pbc\n"
"     The default grid size is 64. Appending a integer value to the\n"
"     command line selects a different value for the grid size.\n"
"     mpirun -np 8 ./poisson_pbc 128 would use a (128,128,128)\n"
"     grid\n\n";

// Required by the BaseThreads class but does not do anything in this program
void *run_threads(void *v) {
}



int main(int argc, char **argv)
{
    int provided, my_rank, npes;
    int nthreads = 1, flag = 0;
    double celldim[6], a0[3], a1[3], a2[3], omega;

    // nthreads won't do anything in this example program but should set it to 1 so as to not waste resources
    BaseThread *B = BaseThread::getBaseThread(nthreads);

    // Set physical grid dimensions but let user override on command line
    int GRIDX = 64, GRIDY = 64, GRIDZ = 64;
    if(argc) {
        if(argv[1]) {
            GRIDX = atoi(argv[1]);
        }
        if(GRIDX < 32) GRIDX = 32;
    }
    GRIDY = GRIDX;
    GRIDZ = GRIDX;

    // MPI node layout needs to be fixed because of the way we are laying out the test charge
    int NODES_X = 2, NODES_Y = 2, NODES_Z = 2;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

    // get my mpi rank and size
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    // get total mpi process count
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    // Check that MPI grid matches npes
    if(npes != (NODES_X * NODES_Y * NODES_Z))
        rmg_error_handler(__FILE__, __LINE__, "Check the number of MPI procs you have requested and make sure it matches (NODES_X * NODES_Y * NODES_Z)");


    // Print Header
    if(my_rank == 0) cout << header;

    // Instantiate and initialize a grid object with a default fine/coarse ratio of 1
    BaseGrid *G = new BaseGrid(GRIDX, GRIDY, GRIDZ, NODES_X, NODES_Y, NODES_Z, 1);

    // Finish setting up the grid object
    G->set_nodes(my_rank);    

    int pbasis = G->get_P0_BASIS(1);

    // Set up a simple cubic cell 20.0 AU per side
    Lattice L;
    L.set_ibrav_type(CUBIC_PRIMITIVE);
    celldim[0] = 20.0;
    celldim[1] = 20.0;
    celldim[2] = 20.0;
    celldim[3] = 0.0;
    celldim[4] = 0.0;
    celldim[5] = 0.0;
    L.latgen(celldim, a0, a1, a2, &omega, &flag);



    // Setup for grids and threads must be completed before we create a TradeImages object
    TradeImages T(G);
    T.set_MPI_comm(MPI_COMM_WORLD);

    // Allocate some arrays to hold the charge and the potential
    double *rho = new double[pbasis];
    double *vh = new double[pbasis];
   
 
    // Initialize them to zero on all nodes
    for(int idx = 0;idx < pbasis;idx++) rho[idx] = 0.0;

    // We need some charge density to work with. For periodic boundary conditions the total cell
    // must be charge neutral otherwise Poissons equation diverges. For this example we use Gaussians
    // centered in each domain containing positive charges simulating ionic cores. We then add a
    // uniform background charge density that simulates the electrons and neutralizes the overall cell.
    int ix, iy, iz, incx, incy;
    double beta=2.0, factor=1.0, rsquared, grid_spacing, vel, totalcharge = 0.0;
    incx = G->get_PY0_GRID(1) * G->get_PZ0_GRID(1);
    incy = G->get_PZ0_GRID(1);
    int icx = G->get_PX0_GRID(1)/2;
    int icy = G->get_PY0_GRID(1)/2;
    int icz = G->get_PZ0_GRID(1)/2;
    grid_spacing = G->get_hxgrid(1) * L.get_xside();      // Cubic so isotropic in all directions
    vel = L.get_omega() / G->get_GLOBAL_BASIS(1);

    for(ix = 0;ix < G->get_PX0_GRID(1);ix++) {
        for(iy = 0;iy < G->get_PY0_GRID(1);iy++) {
            for(iz = 0;iz < G->get_PZ0_GRID(1);iz++) {

                rsquared = grid_spacing * grid_spacing * double((ix-icx)*(ix-icx) + (iy-icy)*(iy-icy) + (iz-icz)*(iz-icz));
                rho[ix*incx + iy*incy + iz] = factor * exp(-beta*rsquared);
                totalcharge += rho[ix*incx + iy*incy + iz] * vel;

            }

        }
    }

    // now neutralize it with a uniform charge
    totalcharge = RmgSumAll(totalcharge, T.get_MPI_comm());
    cout << "For Rank " << my_rank << " total charge = " << totalcharge << endl;
    double t1 = totalcharge / ((double) (npes * pbasis)) / vel;
    for(ix=0;ix < pbasis;ix++) rho[ix] -= t1;

    // Check it. Sum over all MPI nodes should be zero.
    totalcharge = 0.0;
    for(ix=0;ix < pbasis;ix++) totalcharge+=rho[ix] * vel;

    totalcharge = RmgSumAll(totalcharge, T.get_MPI_comm());
    if(my_rank == 0)
        cout << "Total charge for full cell = " << totalcharge << endl;


    // Wait until everyone gets here
    MPI_Barrier(MPI_COMM_WORLD);

    // Set some reasonable multigrid parameters
    int min_sweeps = 5;
    int max_sweeps = 25;
    int global_presweeps = 3;
    int global_postsweeps = 3;
    int mucycles = 2;
    double rms_target = 1.0e-08;
    double global_step = 0.6;
    double coarse_step = 0.6;
    int boundaryflag = PERIODIC;

    // Show convergence acceleration from adding levels
    if(my_rank == 0) cout << endl;

    for(int levels = 0;levels < 5;levels++) {

        double residual;

        // Reinitialize hartree potential to zero before entering solver.
        for(int idx = 0;idx < G->get_PX0_GRID(1) * G->get_PY0_GRID(1) * G->get_PZ0_GRID(1);idx++) vh[idx] = 0.0;

        if(my_rank == 0)
            cout << "Solving poissons equation with " << levels << " multigrid levels. Max " << max_sweeps << " sweeps or RMS target < " << rms_target << endl << "  ";


        residual = CPP_get_vh (G, &L, &T, rho, vh, min_sweeps, max_sweeps, levels, global_presweeps,
                global_postsweeps, mucycles, rms_target,
                global_step, coarse_step, boundaryflag, 1, true);

        if(my_rank == 0)
            cout << endl;
    }

    if(my_rank == 0) cout << endl;

    // Shutdown
    MPI_Finalize ();
}
