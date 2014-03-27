/****f* QMD-MGDFT/get_vh.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void get_vh(double *rho, double *rhoc, double *vh_eig, int sweeps, int maxlevel)
 *   Iterative procedure to obtain the hartree potential.
 *   Uses Mehrstallen finite differencing with multigrid accelerations.
 *   The multigrid scheme uses a standard W-cycle with Jacobi relaxation.
 *   (Underrelaxed on the coarser grids.)
 * INPUTS
 *   rho: total charge density
 *   rhoc: compensating charge density
 *   sweeps:  number of iterations 
 *   maxlevel:  maximum multigrid level
 * OUTPUT
 *   vh_eig: Hartree potential
 * PARENTS
 *   init.c scf.c
 * CHILDREN
 *   getpoi_bc.c pack_vhstod.c pack_ptos.c app_cir.c 
 *   set_bc.c app_cil.c mgrid_solv.c pack_vhdtos.c
 * SOURCE
 */



#include <iostream>

#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "FineGrid.h"
#include "RmgSumAll.h"
#include "auxiliary.h"
#include "vhartree.h"
#include "rmg_error.h"

using namespace std;

#define MAX_MG_LEVELS 8


/// Poisson solver that uses compact implicit (Mehrstellen) and multigrid techniques.
/// @param rho Charge density. When using periodic boundary conditions the cell must be charge neutral.
/// @param vhartree Hartree potential corresponding to rho.
/// @param min_sweeps Minimum number of top level sweeps to perform.
/// @param max_sweeps Maximum number of top level sweeps to perform.
/// @param max_level Maximum number of multigrid levels to use.
/// @param global_presweeps Number of presweeps to use on the finest (0th) grid level.
/// @param global_postsweeps Number of postsweeps to use on the finest (0th) grid level.
/// @param mucycles Number of mu cycles (also known as W-cycles) to use in the multigrid solver.
/// @param rms_target Value for the root mean square residual at which to stop.
/// @param global_step Time step for the jacobi iteration on the finest (0th) grid level.
/// @param coarse_step Time step for the jacobi iteration on the coarse grid levels.
/// @param boundaryflag Type of boundary condition. Periodic is implemented internally.
double CPP_get_vh (double * rho, double *vhartree,
                 int min_sweeps, int max_sweeps, int maxlevel, 
                 int global_presweeps, int global_postsweeps, int mucycles, 
                 double rms_target, double global_step, double coarse_step, int boundaryflag)
{

    int idx, its, cycles;
    double t1, vavgcor, diag;
    double *mgrhsarr, *mglhsarr, *mgresarr, *work;
    double *sg_rho, *sg_vh, *sg_res, *nrho,  residual = 100.0;
    double k_vh;
    FineGrid FG(2);     // Hartree potential is on a double density grid
    Lattice L;
    Mgrid MG(&L);
    TradeImages T;
    int global_basis = FG.get_GLOBAL_BASIS();

    /* Pre and post smoothings on each level */
    int poi_pre[MAX_MG_LEVELS] = { 0, 3, 3, 3, 3, 3, 3, 3 };
    int poi_post[MAX_MG_LEVELS] = { 0, 3, 3, 3, 3, 3, 3, 3 };

    if(maxlevel >= MAX_MG_LEVELS)
       rmg_error_handler(__FILE__, __LINE__, "Too many multigrid levels requested.");

    int dimx = FG.get_PE_GRIDX(), dimy = FG.get_PE_GRIDY(), dimz = FG.get_PE_GRIDZ();

    // Solve to a high degree of precision on the coarsest level
    poi_pre[maxlevel] = 25;

    // Multigrid solver can handle both poisson and helmholtz set k_vh=0 for poisson
    k_vh = 0.0;

    int nits = global_presweeps + global_postsweeps;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);


    /* Grab some memory for our multigrid structures */
    mgrhsarr = new double[sbasis];
    mglhsarr = new  double[sbasis];
    mgresarr = new  double[sbasis];
    work = new  double[4*sbasis];
    sg_rho = new  double[sbasis];
    sg_vh = new  double[sbasis];
    sg_res = new  double[sbasis];
    nrho = new  double[sbasis];

    CPP_app_cir_driver<double> (rho, mgrhsarr, dimx, dimy, dimz, APP_CI_FOURTH);



    /* Multiply through by 4PI */
    t1 = -FOUR * PI;
    for(idx = 0;idx < pbasis;idx++) mgrhsarr[idx] = t1 * mgrhsarr[idx];

    its = 0;

    while ( ((its < max_sweeps) && (residual > rms_target))  || (its < min_sweeps))
    {


        /* Mehrstallen smoothings */
        for (cycles = 0; cycles < nits; cycles++)
        {


	    /*At the end of this force loop, laplacian operator is reapplied to evalute the residual. Therefore,
	     * there is no need to reapply it, when this loop is called second, third, etc time. */
            if ( (cycles) || (!its))
            {

                /* Apply operator */
                diag = CPP_app_cil_driver (vhartree, mglhsarr, dimx, dimy, dimz,
                            FG.get_hxxgrid(), FG.get_hyygrid(), FG.get_hzzgrid(), APP_CI_FOURTH);
                diag = -1.0 / diag;

                /* Generate residual vector */
                for (idx = 0; idx < pbasis; idx++)
                {

                    mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];

                }                   /* end for */

            }

            /* Pre and Post smoothings and multigrid */
            if (cycles == global_presweeps)
            {


                /* Transfer res into smoothing grid */
                CPP_pack_ptos<double> (sg_res, mgresarr, dimx, dimy, dimz);

                MG.mgrid_solv<double> (mglhsarr, sg_res, work,
                            dimx, dimy, dimz, 
                            FG.get_hxxgrid(), FG.get_hyygrid(), FG.get_hzzgrid(),
                            0, FG.get_neighbors(), maxlevel, poi_pre,
                            poi_post, mucycles, coarse_step, k_vh,
                            FG.get_GLOBAL_GRIDX(), FG.get_GLOBAL_GRIDY(), FG.get_GLOBAL_GRIDZ(),
                            FG.get_PE_OFFSETX(), FG.get_PE_OFFSETY(), FG.get_PE_OFFSETZ(),
                            FG.get_PE_GRIDX(), FG.get_PE_GRIDY(), FG.get_PE_GRIDZ(), boundaryflag);


                /* Transfer solution back to mgresarr array */
                CPP_pack_stop<double> (mglhsarr, mgresarr, dimx, dimy, dimz);

                /* Update vh */
                t1 = ONE;
                for(int i = 0;i < pbasis;i++) vhartree[i] += mgresarr[i];

            }
            else
            {

                /* Update vh */
                t1 = - global_step * diag;
                for(int i = 0;i < pbasis;i++) vhartree[i] = vhartree[i] + t1 * mgresarr[i];

            }                   /* end if */

            if (boundaryflag == PERIODIC)
            {

                /* Evaluate the average potential */
                vavgcor = 0.0;
                for (idx = 0; idx < pbasis; idx++)
                    vavgcor += vhartree[idx];

                vavgcor =  RmgSumAll(vavgcor, T.get_MPI_comm());
                t1 = (double) global_basis;
                vavgcor = vavgcor / t1;


                /* Make sure that the average value is zero */
                for (idx = 0; idx < pbasis; idx++)
                    vhartree[idx] -= vavgcor;

            }                   /* end if */

        }                       /* end for */

        /*Get residual*/
        diag = CPP_app_cil_driver<double> (vhartree, mglhsarr, dimx, dimy, dimz,
                            FG.get_hxxgrid(), FG.get_hyygrid(), FG.get_hzzgrid(), APP_CI_FOURTH);
        diag = -1.0 / diag;
        residual = 0.0;

        /* Generate residual vector */
        for (idx = 0; idx < pbasis; idx++)
        {

            mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];
            residual += mgresarr[idx] * mgresarr[idx];

        }                   /* end for */

        residual = sqrt (RmgSumAll(residual, T.get_MPI_comm()) / global_basis);

        //cout << "\n get_vh sweep " << its << " rms residual is " << residual;

	    
        its ++;
    }                           /* end for */

//    if(FG.get_gridpe() == 0)
//        cout << "get_vh: executed " << its << " sweeps, residual is " << residual << endl;

    /* Release our memory */
    delete [] nrho;
    delete [] sg_res;
    delete [] sg_vh;
    delete [] sg_rho;
    delete [] work;
    delete [] mgresarr;
    delete [] mglhsarr;
    delete [] mgrhsarr;

    return residual;

} // end CPP_get_vh


