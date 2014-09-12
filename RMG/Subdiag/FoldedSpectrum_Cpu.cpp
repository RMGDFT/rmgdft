#include <complex>
#include <omp.h>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

// I have not finished updating this to work with complex orbitals yet. Given that the folded spectrum method is only
// useful for large systems which are almost always run at gamma with real orbitals it's not a high priority but should
// be straightforward enough to finish.
template int FoldedSpectrumCpu<double> (Kpoint<double> *, int, double *, int, double *, int, double *, double *, int, int *, int, double *);

template <typename KpointType>
int FoldedSpectrumCpu(Kpoint<KpointType> *kptr, int n, KpointType *A, int lda, KpointType *B, int ldb, 
		double *eigs, double *work, int lwork, int *iwork, int liwork, KpointType *C)
{

    RmgTimer RT0("Diagonalization: fs:");

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    char *trans="n", *transt="t", *transn="n";
    char *cuplo = "l", *side="l", *diag="n", *jobz="V";

    int ione=1, itype=1, info=0;
    double rone = 1.0;
    double alpha, beta=0.0, lambda;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;


    int NPES = G->get_PE_X() * G->get_PE_Y() * G->get_PE_Z();

    // Array storage for folded spectrum diagonalization
    static int *fs_eigstart = NULL;
    static int *fs_eigstop = NULL;
    static int *fs_eigcounts = NULL;

    // Initialize them
    if(!fs_eigstart) {
        fs_eigstart = new int[NPES];
        fs_eigstop = new int[NPES];
        fs_eigcounts = new int[NPES];
    }

    for(int idx = 0;idx < NPES;idx++) {
        double t1 = (double)n;
        t1 = t1 / ((double)NPES);
        double t2 = t1 * (double)idx;
        fs_eigstart[idx] = (int)rint(t2);
        fs_eigstop[idx] = (int)rint(t1 + t2);
        fs_eigcounts[idx] = fs_eigstop[idx] - fs_eigstart[idx];
        fs_eigstart[idx] *= n;
        fs_eigstop[idx] *= n;
        fs_eigcounts[idx] *= n;

    }

     
    // Folded spectrum method is parallelized over PE's. Each PE gets assigned
    // a subset of the eigenvectors.
    double t1 = (double)n;
    t1 = t1 / ((double)NPES);
    double t2 = t1 * (double)pct.gridpe;
    int eig_start = (int)rint(t2);
    int eig_stop = (int)rint(t1 + t2);
    int eig_step = eig_stop - eig_start;
    if(pct.gridpe == (NPES - 1)) eig_stop = n;

    // Set width of window in terms of a percentage of n. Larger values will be slower but
    // exhibit behavior closer to full diagonalization.
    double r_width = 0.3;
    t1 = (double)n;
    int n_win = (int)(r_width * t1);

    // Find start of interval
    int ix = n_win - eig_step;
    if(ix < 4)
        rmg_error_handler(__FILE__, __LINE__, "Too few PE's to use folded spectrum method for this problem");
    if(ix % 2) {
        n_win++;
        ix = n_win - eig_step;
    }

    int n_start = eig_start - ix/2;
    if(n_start < 0)n_start = 0;
    if((n_start + n_win) > n) {
        n_start = n - n_win;
    }

    double *Vdiag = new double[n];
    double *p0 = new double[n];
    double *p1 = new double[n];
    double *d_p0 = new double[n * eig_step];
    double *d_p1 = new double[n];
    double *tarr = new double[n];
    double *Asave = new double[n*n];



    RmgTimer *RT1 = new RmgTimer("Diagonalization: fs: cholesky");
    //  Form a Cholesky factorization of B.
    dpotrf(cuplo, &n, B, &ldb, &info);
    if( info != 0 )
        rmg_error_handler(__FILE__, __LINE__, "dpotrf failure");
    delete(RT1);


    //  Transform problem to standard eigenvalue problem
    RT1 = new RmgTimer("Diagonalization: fs: tridiagonal");
    dsygst_(&itype, cuplo, &n, A, &lda, B, &ldb, &info);
    if( info != 0 ) 
        rmg_error_handler(__FILE__, __LINE__, "dsygst failure");
    delete(RT1);


    // We need to wait until a is diagonally dominant so we skip the first 3 steps
    if((ct.scf_steps < 3) && (ct.runflag != RESTART)) {

        RT1 = new RmgTimer("Diagonalization: fs: init");
        dsyevd(jobz, cuplo, &n, A, &lda, eigs,
                         work,  &lwork,
                         iwork, &liwork,
                         &info);
        if( info != 0 ) 
            rmg_error_handler(__FILE__, __LINE__, "dsyevd failure");

        //   For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        //        backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
        //   dtrsm_( "Leftx", uplo, trans, "Non-unit", &n, &n, &rone, b, &ldb, a, &lda );
        //
        dtrsm (side, cuplo, transt, diag, &n, &n, &rone, B, &ldb, A, &lda);
        delete(RT1);

    }
    else {
       
        RT1 = new RmgTimer("Diagonalization: fs: normal");
        // Zero out matrix of eigenvectors (V) and eigenvalues n. G is submatrix storage
        KpointType *V = new KpointType[n*n]();
        KpointType *G = new KpointType[n*n]();
        double *n_eigs = new double[n]();

        // AX=lambdaX  store a copy of A in Asave
        for(int idx = 0;idx < n*n;idx++) Asave[idx] = A[idx];
        //QMD_dcopy (n*n, a, 1, Asave, 1);

     
        // Do the submatrix along the diagonal to get starting values for folded spectrum
        //--------------------------------------------------------------------
        RmgTimer *RT2 = new RmgTimer("Diagonalization: fs: submatrix");
        for(int ix = 0;ix < n_win;ix++){
            for(int iy = 0;iy < n_win;iy++){
                G[ix*n_win + iy] = Asave[(n_start+ix)*n + n_start + iy];
            }
        }

        for(int idx = 0;idx < n_win * n_win;idx++) A[idx] = G[idx];
        //QMD_dcopy (n_win * n_win, G, 1, a, 1);

        dsyevd_(jobz, cuplo, &n_win, A, &n_win, &eigs[n_start],
                        work, &lwork,
                        iwork, &liwork,
                        &info);
        if( info != 0 ) 
                rmg_error_handler(__FILE__, __LINE__, "dsyevd failure");

        //--------------------------------------------------------------------

        for(int idx = 0;idx < n_win * n_win;idx++) G[idx] = A[idx];
        //QMD_dcopy (n_win * n_win, a, 1, G, 1);

        for(int ix = 0;ix < n_win;ix++) {
            Vdiag[ix] = ONE_t;
            if(G[ix*n_win + ix] < 0.0) Vdiag[ix] = -ONE_t;
        }

        // Store the eigen vector from the submatrix
        for(int ix = 0;ix < n_win;ix++) {

            if(((n_start+ix) >= eig_start) && ((n_start+ix) < eig_stop)) {

                for(int iy = 0;iy < n_win;iy++) {
                      V[(ix + n_start)*n + n_start + iy] = Vdiag[ix] * G[ix * n_win + iy];
                }

            }

        }
        delete(RT2);


        // Apply folded spectrum to this PE's range of eigenvectors
        RT2 = new RmgTimer("Diagonalization: fs: iteration");
#if 0
        for(int eig_index = eig_start;eig_index < eig_stop;eig_index++) {

                lambda = eigs[eig_index];
                n_eigs[eig_index] = lambda;

                for(int ix=0;ix<n;ix++){
                    Asave[ix*n + ix] -= lambda;
                }
                alpha = 1.0;
                for(int idx = 0;idx < n*n;idx++) C[idx] = Asave[idx];
                //QMD_dcopy (n*ct.num_states, Asave, 1, C, 1);

                // Restore matrix for next pass
                for(ix=0;ix<n;ix++){
                    Asave[ix*n + ix] += lambda;
                }

                for(int idx = 0;idx < n;idx++) d_p0[idx] = V[eig_index*n + idx];
                //QMD_dcopy (n, &V[eig_index*n], 1, d_p0, 1);

                alpha = -0.5;
                beta = 0.0;
                for(int its = 0;its < 6;its++) {
                    dgemv_(transn, &n, &n, &alpha, C, &n, d_p0, &ione, &beta, d_p1, &ione);
                    daxpy(&n, &rone, d_p1, &ione, d_p0, &ione);
                }
#endif
        for(int eig_index = eig_start;eig_index < eig_stop;eig_index++) {
            n_eigs[eig_index] = eigs[eig_index];
            int offset = n * (eig_index - eig_start);
            for(int idx = 0;idx < n;idx++) d_p0[offset + idx] = V[eig_index*n + idx];
        }

        FoldedSpectrumIterator(Asave, n, &eigs[eig_start], eig_stop - eig_start, d_p0, -0.5, 6);


        for(int eig_index = eig_start;eig_index < eig_stop;eig_index++) {

                int offset = n * (eig_index - eig_start);
                // Renormalize
                //t1 = dnrm2_(&n, d_p0, &ione);
                t1 = 0.0;
                for(int idx = 0;idx < n;idx++) t1 += std::norm(d_p0[offset + idx]);
                t1 = 1.0 / sqrt(t1);
                for(int idx = 0;idx < n;idx++) d_p0[offset + idx] *= t1;
                //dscal(&n, &t1, d_p0, &ione);

                for(int idx = 0;idx < n;idx++) V[eig_index*n + idx] = d_p0[offset + idx];
                //QMD_dcopy (n, d_p0, 1, &V[eig_index*n], 1);
        }
        delete(RT2);


        // Make sure all PE's have all eigenvectors. Possible optimization here would be to 
        // overlap computation in the above loop with communication here.
//        MPI_Allreduce(MPI_IN_PLACE, V, n*n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        RT2 = new RmgTimer("Diagonalization: fs: allreduce1");
        MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);
        delete(RT2);

        // Do the same for the eigenvalues
        // Could replace this with an MPI_Allgatherv as well
        RT2 = new RmgTimer("Diagonalization: fs: allreduce2");
        MPI_Allreduce(MPI_IN_PLACE, n_eigs, n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        delete(RT2);

        // Copy summed eigs back to w
        for(int idx = 0;idx < n;idx++) eigs[idx] = n_eigs[idx];

        // Gram-Schmidt ortho for eigenvectors.
        RT2 = new RmgTimer("Diagonalization: fs: Gram-Schmidt");
        alpha = 1.0;
        beta = 0.0;

        // Overlaps
        RmgTimer *RT3 = new RmgTimer("Diagonalization: fs: overlaps");
        //QMD_dcopy (n*n, V, 1, a, 1);
        for(int idx = 0;idx < n*n;idx++) A[idx] = V[idx];
        dsyrk (cuplo, transt, &n, &n, &alpha, A, &n, &beta, C, &n);
        delete(RT3);

        // Cholesky factorization
        RT3 = new RmgTimer("Diagonalization: fs: cholesky");
        dpotrf(cuplo, &n, C, &n, &info);
        delete(RT3);

        // Get inverse of diagonal elements
        for(ix = 0;ix < n;ix++) tarr[ix] = 1.0 / C[n * ix + ix];

//----------------------------------------------------------------
        for(int idx = 0;idx < n*n;idx++)G[idx] = ZERO_t;

        int idx, omp_tid;
        double *darr, *sarr;
        int st, st1;
#pragma omp parallel private(idx,st,st1,omp_tid,sarr)
{
        omp_tid = omp_get_thread_num();
        if(omp_tid == 0) darr = new double[n * omp_get_num_threads()];
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
        for(idx = eig_start;idx < eig_stop;idx++) {

            sarr = &darr[omp_tid*n];

            for (int st = 0; st < n; st++) sarr[st] = V[st*n + idx];

            for (int st = 0; st < n; st++) {

                sarr[st] *= tarr[st];

                for (int st1 = st+1; st1 < n; st1++) {
                    sarr[st1] -= C[st1 + n*st] * sarr[st];
                }

            }

            for (st = 0; st < n; st++) G[st*n + idx] = sarr[st];

        }
} // end omp section

        delete [] darr;

        delete(RT2);


        // A matrix transpose here would let us use an Allgatherv which would be
        // almost twice as fast for the communications part.
        RT2 = new RmgTimer("Diagonalization: fs: allreduce3");
        MPI_Allreduce(MPI_IN_PLACE, G, n*n * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        delete(RT2);
        for(int idx = 0;idx < n*n;idx++) A[idx] = G[idx];
        //QMD_dcopy (n*n, G, 1, a, 1);



        RT2 = new RmgTimer("Diagonalization: fs: dtrsm");
        dtrsm (side, cuplo, transt, diag, &n, &n, &rone, B, &ldb, A, &lda);
        delete(RT2);

        delete(RT1);
        delete [] V;
        delete [] G;
        delete [] n_eigs;

    }


    delete [] Asave;
    delete [] tarr;
    delete [] d_p1;
    delete [] d_p0;
    delete [] p1;
    delete [] p0;
    delete [] Vdiag;

    return 0;
} 
