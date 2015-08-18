/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
void DiagElemental(int n, double *H, double *S)
{

    // MPI_COMM_WORLD. There is another constructor that allows you to 
    // specify the grid dimensions, Grid g( comm, r ), which creates a
    // grid of height r.
    Grid g( mpi::COMM_WORLD );

    printf("\n set Comm_g");
    fflush(NULL);
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
    int tem[n+npes];
    int num_local = (n + npes -1)/npes;
    for (int i = 0; i < npes; i++)
    {
        for(int j = 0; j < num_local;j++)
            tem[i + j * npes] = i * num_local + j;
    }
    
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        // Our process owns the rows colShift:colStride:n,
        //           and the columns rowShift:rowStride:n
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
     //       printf("\n aaaa  %d %d %d %d %f", i, j, iLoc, jLoc, H[iLoc*n+tem[jLoc]]);
            A.SetLocal( iLoc, jLoc, H[iLoc*n+tem[jLoc]]);
            B.SetLocal( iLoc, jLoc, S[iLoc*n+tem[jLoc]]);
        }
    }

    // Call the eigensolver. We first create an empty complex eigenvector 
    // matrix, X[MC,MR], and an eigenvalue column vector, w[VR,* ]
    //
    // Optional: set blocksizes and algorithmic choices here. See the 
    //           'Tuning' section of the README for details.
    DistMatrix<double,VR,STAR> w(g );
    DistMatrix<double, VC, STAR> X(n,n, g );
    //DistMatrix<double> X(g);

    HermitianEigSubset<double> subset;
    HermitianEigCtrl<double> ctrl_d;
       // ctrl_d.timeStages = true;
       // ctrl_d.tridiagCtrl.symvCtrl.bsize = 128;
       // ctrl_d.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = true;
    //ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;

      ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_d.tridiagCtrl.order = COLUMN_MAJOR;

    //HermitianGenDefEig(AXBX, LOWER, A, B, w, X, ASCENDING, subset, ctrl_d); 
    SetBlocksize(128);
    HermitianGenDefEig(AXBX, LOWER, A, B, w, X, ASCENDING, subset, ctrl_d);

    Print(w, "w");
 //   for(int i = 0; i < n; i++)
 //   {
 //       printf(" %f  ", w.GetRealPart(i,i) * 27.2);
 //       if(i%5==0) printf("\n");
 //   }
    fflush(NULL);
    //HermitianEig(LOWER, A, w, X, ASCENDING ); 


}
