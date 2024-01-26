/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/


#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"
#include "common_prototypes1.h"

#if ELEMENTAL_LIBS


#include "El.hpp"
using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
void DiagElemental(STATE *states, int n, double *H, double *S, 
        double *rho_matrix, double *theta_ptr)
{

    // MPI_COMM_WORLD. There is another constructor that allows you to 
    // specify the grid dimensions, Grid g( comm, r ), which creates a
    // grid of height r.
    RmgTimer *RT= new RmgTimer("3-DiagElemental");
    Grid g( mpi::COMM_WORLD );

    // Create an n x n complex distributed matrix, 
    // We distribute the matrix using grid 'g'.
    //
    // There are quite a few available constructors, including ones that 
    // allow you to pass in your own local buffer and to specify the 
    // distribution alignments (i.e., which process row and column owns the
    // top-left element)
    DistMatrix<double, VC, STAR> A( n, n, g );
    DistMatrix<double, VC, STAR> B( n, n, g );
    // DistMatrix<double> A( n, n, g );
    // DistMatrix<double> B( n, n, g );

    // Fill the matrix since we did not pass in a buffer. 
    //
    // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
    // the global matrix is Hermitian. However, only one triangle of the 
    // matrix actually needs to be filled, the symmetry can be implicit.
    //
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    int npes = mpi::Size(mpi::COMM_WORLD);
    int *perm;
    int i,j, index_restart, index_pe;
    int num_local = (n + npes -1)/npes;
        


    perm = new int[n+npes];
    index_restart =0;
    for (i = 0; i < npes; i++)
    {
        for(j = 0; j < num_local;j++)
        {

            if (j * npes + i >= n)
            {
                index_restart = i*num_local +j;
                index_pe = i+1;
                i+= npes;
                break;
            }
            perm[j * npes + i]=  i * num_local + j;
        }
    }
    if(index_restart >0)
    {

        for (i = index_pe; i < npes; i++)
        {
            for(j = 0; j < num_local-1; j++)
            {
                int item= index_restart + (i-index_pe) * (num_local-1) +j;
                perm[j * npes + i] = item;
            }
        }
    }


    //    if(g.Rank() == 0) 
    //        for(int i = 0; i < n; i++) rmg_printf("\n aaaxx %d  %d\n", i, perm[i]);

    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            A.SetLocal( iLoc, jLoc, H[iLoc*n+perm[jLoc]]);
            B.SetLocal( iLoc, jLoc, S[iLoc*n+perm[jLoc]]);
        }
    }

    // Optional: set blocksizes and algorithmic choices here. See the 
    //           'Tuning' section of the README for details.
    DistMatrix<double,VR,STAR> w(g );
    DistMatrix<double, VC, STAR> X(n,n, g );

    HermitianEigSubset<double> subset;
    HermitianEigCtrl<double> ctrl_d;
    // ctrl_d.timeStages = true;
    // ctrl_d.tridiagCtrl.symvCtrl.bsize = 128;
    // ctrl_d.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = true;
    //ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;

    ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl_d.tridiagCtrl.order = COLUMN_MAJOR;
    ctrl_d.useSDC= true;   //  does not work
    ctrl_d.sdcCtrl.tol = 1.0e-15;

    SetBlocksize(128);
    RmgTimer *RT1= new RmgTimer("3-DiagElemental: AXBX");
   // HermitianGenDefEig(AXBX, LOWER, A, B, w, X, ASCENDING, subset, ctrl_d); 
    HermitianGenDefEig(AXBX, LOWER, A, B, w, X, ASCENDING);
    delete(RT1);

    //Print(w, "w");

    for (int st1 = 0; st1 < n; st1++)
    {
        states[st1].eig[0] = w.Get(st1, 0);
    }



    //overwrite B with X[i,j] * occ[j]


    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            B.SetLocal( iLoc, jLoc, X.GetLocal(iLoc, jLoc) * states[j].occupation[0]);
        }
    }


    RmgTimer *RT2= new RmgTimer("3-DiagElemental: Gemm");
    Gemm( NORMAL, TRANSPOSE, 1., B, X, 0., A);
    delete(RT2);

    // A = X *occ * X

    
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            rho_matrix[iLoc*n+perm[jLoc]]= A.GetLocal( iLoc, jLoc);
        }
    }
    delete(RT);

    RmgTimer *RT3= new RmgTimer("3-mg_eig: (S^-1)H");
    //  calculating  S^-1 H
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            A.SetLocal( iLoc, jLoc, H[iLoc*n+perm[jLoc]]);
            B.SetLocal( iLoc, jLoc, S[iLoc*n+perm[jLoc]]);
        }
    }

// HPD Solve is a little faster than SymmetricSolve
    //SymmetricSolve(LOWER, NORMAL, B,A);
    HPDSolve(LOWER, NORMAL, B,A);

    Transpose(A, B);
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            theta_ptr[iLoc*n+perm[jLoc]]= 2.0 * B.GetLocal( iLoc, jLoc);
        }
    }
    delete(RT3);

    A.Empty();
    B.Empty();
    X.Empty();
    w.Empty();

}
#endif
