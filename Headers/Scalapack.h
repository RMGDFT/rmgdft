/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#ifndef RMG_Scalapack_H
#define RMG_Scalapack_H 1


#include <mpi.h>
#include <complex>
#include "rmg_mangling.h"


/* Blacs dimension */
#define DLEN    9

// Maximum number of scalapack groups
#define MAX_SCALAPACK_GROUPS 32
// Maximum number of folded spectrum scalapack instances
#define MAX_FOLDED_SCALAPACKS 32


#ifdef __cplusplus

// npes is the total number of pes, ngroups is the number of
// scalapack groups that should be created out of npes where
// npes <= group_pes * ngroups * images_per_node

class Scalapack {

public:

    explicit Scalapack(int ngroups, int thisimg, int images_per_node, int N, int NB, int last, MPI_Comm rootcomm);
    void DistributeMatrix(double *A, double *A_dist);
    void GatherMatrix(double *A, double *A_dist);

    void DistributeMatrix(std::complex<double> *A, std::complex<double> *A_dist);
    void GatherMatrix(std::complex<double> *A, std::complex<double> *A_dist);

    int GetRootNpes(void);
    int GetNumGroups(void);
    int GetNumGroupsNext(void);
    int GetScalapackNpes(void);
    int GetGroupIndex(void);
    int GetRows(void);
    int GetCols(void);
    int GetRow(void);
    int GetCol(void);
    int GetDistMdim(void);
    int GetDistNdim(void);
    int ComputeMdim(int m);
    int ComputeNdim(int n);
    int GetCommRank(void);
    int GetRootRank(void);
    int *GetDistDesca(void);
    int ComputeDesca(int m, int n, int *desca);
    int GetIpivSize(void);
    int GetContext(void);
    bool Participates(void);
    Scalapack *GetNextScalapack(void);
    MPI_Comm GetComm(void);
    MPI_Comm GetRootComm(void);

    void Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op);
    void ScalapackBlockAllreduce(double *buf, size_t count);
    void BcastRoot(void *buffer, int count, MPI_Datatype datatype);
    void BcastComm(void *buffer, int count, MPI_Datatype datatype);

    void Pgemm (char *transa, char *transb, int *M, int *N, int *K, double *alpha,
                       double *A, int *IA, int *JA, int *desca,
                       double *B, int *IB, int *JB, int *descb,
                       double *beta, double *C, int *IC, int *JC, int *descc);

    void Pgemm (char *transa, char *transb, int *M, int *N, int *K, std::complex<double> *alpha,
                       std::complex<double> *A, int *IA, int *JA, int*desca,
                       std::complex<double> *B, int *IB, int *JB, int *descb,
                       std::complex<double> *beta, std::complex<double> *C, int *IC, int *JC, int *descc);

    void Pgesv (int *N, int *NRHS, double *A, int *IA, int *JA, int *desca, int *ipiv, double *B, int *IB,
                            int *JB, int *descb, int *info);

    void Pgesv (int *N, int *NRHS, std::complex<double> *A, int *IA, int *JA, int *desca, int *ipiv, std::complex<double> *B, int *IB,
                            int *JB, int *descb, int *info);

    void CopySquareMatrixToDistArray(double *A, double *A_dist, int n, int *desca);
    void CopySquareMatrixToDistArray(std::complex<double> *A, std::complex<double> *A_dist, int n, int *desca);

    void CopyDistArrayToSquareMatrix(double *A, double *A_dist, int n, int *desca);
    void CopyDistArrayToSquareMatrix(std::complex<double> *A, std::complex<double> *A_dist, int n, int *desca);

    // Next level of scalapack
    Scalapack *next;

    size_t *dist_sizes;    // sizes of distributed matrices on all PEs
    size_t *dist_offsets;  // offsets into global matrix of each PE distributed matrix

    ~Scalapack(void);

protected:

    void matscatter (double *globmat, double *dismat, int size, int *desca, bool isreal);
    void matgather (double *globmat, double *dismat, int size, int *desca, bool isreal);


    int N;              // Operates on matrices of size (N,N)
    int scalapack_npes; // number of processors from the root_com that participate in sp operatrions
    int NB;             // Blocking factors
    int context;        // blacs context of this group of pes
    int npes;           // total number of pes in the root_comm
    int ngroups;        // total number of groups
    int ngroups_next;   // total number of groups at the next level
    int group_pes;      // number of pes in this group
    int group_index;    // index of this group
    int group_rows;     // rows in this group
    int group_cols;     // cols in this group
    int my_row;         // blacs row of this PE
    int my_col;         // blacs col of this PE
    int root_rank;      // rank in rootcomm
    int comm_rank;      // rank in this comm
    int *local_desca;   // desca for local matrix
    int *dist_desca;    // desca for distributed matrix
    bool participates;  // whether or not this PE participates in scalapack calculations
    int m_dist;         // rows of distributed matrix
    int n_dist;         // cols of distributed matrix
    MPI_Comm comm;      // communicator for a specific group
    MPI_Comm root_comm; // parent communicator
    MPI_Comm used_comm; // holds processes that participate in scalapack ops
                        // may include more than one group
    MPI_Comm broadcast_comm;  // holds root process in group 0 and any leftovers

};

#endif


#define		numroc		RMG_FC_GLOBAL(numroc, NUMROC)
#define		indxg2p		RMG_FC_GLOBAL(indxg2p, INDXG2P)
#define		descinit	RMG_FC_GLOBAL(descinit, DESCINIT)
#define		pdsyev		RMG_FC_GLOBAL(pdsyev, PDSYEV)
#define		pcheev		RMG_FC_GLOBAL(pcheev, PCHEEV)
#define		pspocon		RMG_FC_GLOBAL(pspocon, PSPOCON)
#define		pdgemm		RMG_FC_GLOBAL(pdgemm, PDGEMM)
#define		pdlaset		RMG_FC_GLOBAL(pdlaset, PDLASET)
#define		pzgemm		RMG_FC_GLOBAL(pzgemm, PZGEMM)
#define		pdgesv		RMG_FC_GLOBAL(pdgesv, PDGESV)
#define		pzgesv		RMG_FC_GLOBAL(pzgesv, PZGESV)
#define		psubdiag	RMG_FC_GLOBAL(psubdiag, PSUBDIAG)
#define		pdsygvx		RMG_FC_GLOBAL(pdsygvx, PDSYGVX)
#define		pdsyevx		RMG_FC_GLOBAL(pdsyevx, PDSYEVX)
#define		pdpotrf		RMG_FC_GLOBAL(pdpotrf, PDPOTRF)
#define		pdsyngst	RMG_FC_GLOBAL(pdsyngst, PDSYNGST)
#define		pdtrsm		RMG_FC_GLOBAL(pdtrsm, PDTRSM)
#define		pdtran		RMG_FC_GLOBAL(pdtran, PDTRAN)
#define		pztranc		RMG_FC_GLOBAL(pztranc, PZTRANC)
#define		pzhegvx		RMG_FC_GLOBAL(pzhegvx, PZHEGVX)
#define		pdsyevd		RMG_FC_GLOBAL(pdsyevd, PDSYEVD)
#define		pdgeadd		RMG_FC_GLOBAL(pdgeadd, PDGEADD)
#define		pzgeadd		RMG_FC_GLOBAL(pzgeadd, PZGEADD)
#define		pzpotrf		RMG_FC_GLOBAL(pzpotrf, PZPOTRF)
#define		pzhegst		RMG_FC_GLOBAL(pzhegst, PZHEGST)
#define		pzheevd		RMG_FC_GLOBAL(pzheevd, PZHEEVD)
#define		pztrsm		RMG_FC_GLOBAL(pztrsm, PZTRSM)
#define		pdsygst		RMG_FC_GLOBAL(pdsygst, PDSYGST)
#define		pdsyrk		RMG_FC_GLOBAL(pdsyrk, PDSYRK)
#define		pdsymm		RMG_FC_GLOBAL(pdsymm, PDSYMM)
#define		pdgetrf		RMG_FC_GLOBAL(pdgetrf, PDGETRF)
#define		pzgetrf		RMG_FC_GLOBAL(pzgetrf, PZGETRF)
#define		pzgeqpf		RMG_FC_GLOBAL(pzgeqpf, PZGEQPF)
#define		pdgeqpf		RMG_FC_GLOBAL(pdgeqpf, PDGEQPF)
#define		pzgetri		RMG_FC_GLOBAL(pzgetri, PZGETRI)
#define		pdgetri		RMG_FC_GLOBAL(pdgetri, PDGETRI)
#define		pdgetrs		RMG_FC_GLOBAL(pdgetrs, PDGETRS)
#define		pdpocon		RMG_FC_GLOBAL(pdpocon, PDPOCON)
#define		pztranu		RMG_FC_GLOBAL(pztranu, PZTRANU)
#ifdef __cplusplus
extern "C" {
#endif

int Csys2blacs_handle(MPI_Comm SysCtxt );
MPI_Comm Cblacs2sys_handle (int BlacsCtxt);
void Cpdgemr2d(int, int, double*, int, int, int*, double*, int, int, int*, int);
int numroc (int *, int *, int *, int *, int *);
int indxg2p (int *, int *, int *, int *, int *);
void pdgetrf( int *, int *, double *, int *, int *, int *, int *, int * );
void pdgetrs( char *, int *, int *, double *, int *, int *, int *, int *, double *, int *, int *, int *, int *);
void descinit (int[], int *, int *, int *, int *, int *, int *, int *, int *,
               int *);
void pdgesv (int *, int *, double *, int * , int *, int *, int *, double *,
        int *, int *, int *, int *);
void pzgesv (int *, int *, std::complex<double> *, int * , int *, int *, int *, std::complex<double> *,
        int *, int *, int *, int *);
void pdgemm (char *, char *, int *, int *, int *, double *, double *, int *,
             int *, int *, double *, int *, int *, int *, double *, double *,
             int *, int *, int *);
void pzgemm (char *, char *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *,
             int *, int *, std::complex<double> *, int *, int *, int *, std::complex<double> *, std::complex<double> *,
             int *, int *, int *);
void pdsyev (char *, char *, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, int *);
void pdsyevd (char *, char *, int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *, int *, int *, int *, int *);
void pzheevd (char *, char *, int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *, int *, double *, int *, int *, int *, int *);
void pcheev (char *, char *, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, double *, int *, int *);
void pdsygst(int *, char *, int *, double *, int *, int *, int *, double *, int *,
              int *, int *, double *, int *);
void pdsyngst(int *, char *, int *, double *, int *, int *, int *, double *, int *,
              int *, int *, double *, double *, int *, int *);

void pzhegst(int *, char *, int *, double *, int *, int *, int *, double *, int *,
              int *, int *, double *, int *);
void pdpotrf(char *, int*, double*, int*, int*, int*, int*);
void pdpocon(char *, int*, double*, int*, int*, int*, double *, double *, double *, int *, int *, int *, int *);
void pzpotrf(char *, int*, double*, int*, int*, int*, int*);
void pdtrtri(char *, char *, int*, double*, int*, int*, int*, int*);
void pdsyrk( char *, char *, int *, int *, double *, double *, int *, int *, int *,
             double *, double *, int *, int *, int *);
void pdlaset(char *, int *, int *, double *, double *, double *, int *, int *, int *);
void psubdiag (char *, char *, int, double *, int, double *, int *);
void pdsygvx(int *, char*, char*, char*, int*, double *, int*, int*, int*, double*, int*, int*,
       int*, double*, double *, int*, int*, double*, int*, int*, double*, double*, double*, int*,
       int*, int*, double*, int*, int*, int*, int*, int*, double*, int*);
void pzhegvx(int *, char*, char*, char*, int*, double *, int*, int*, int*, double*, int*, int*,
       int*, double*, double *, int*, int*, double*, int*, int*, double*, double*, double*, int*,
       int*, int*, double*, int *, double *, int*, int*, int*, int*, int*, double*, int*);
void pdsyevx(char*, char*, char*, int*, double *, int*, int*, int*, double*, double*, int*,
       int*, double*, int*, int*, double*, double*, double*, int*,
       int*, int*, double*, int*, int*, int*, int*, int*, double*, int*);
void pdgeadd(char *, int *, int *, double *, double *, int *, int *, int *, double *,
       double *, int *, int *, int *);               
void pzgeadd(char *, int *, int *, double *, double *, int *, int *, int *, double *,
       double *, int *, int *, int *);               
void pdtrmm(char *side, char *uplo, char *trans, char *diag, int * m, int *n, double *alpha,
             double * a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);               
void pdtrsm(char *side, char *uplo, char *trans, char *diag, int * m, int *n, double *alpha,
             double * a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);               
void pztrmm(char *side, char *uplo, char *trans, char *diag, int * m, int *n, double *alpha,
             double * a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);               
void pztrsm(char *side, char *uplo, char *trans, char *diag, int * m, int *n, double *alpha,
             double * a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);               
void pdtran( int * M, int * N, double * ALPHA, double * A, int * IA, int * JA, int * DESCA,
              double * BETA, double * C, int * IC, int * JC, int * DESCC );
void pdsymm(char *SIDE, char *UPLO, int *M, int *N, double *ALPHA, double *A, int *IA, int *JA, int *DESCA,
            double *B, int *IB, int *JB, int *DESCB, double *BETA, double *C, int *IC, int *JC, int *DESCC);
void pztranc(int *M, int *N, std::complex<double> *ALPHA, std::complex<double> *A, int *IA, int *JA, int *DESCA, std::complex<double> *BETA, std::complex<double> *C, int *IC, int *JC, int *DESCC);
void pztranu (int *M, int *N, std::complex<double> *ALPHA, std::complex<double> *A, int *IA, int *JA, int *DESCA, std::complex<double> *BETA, std::complex<double> *C, int *IC, int *JC, int *DESCC);
void pzgetri(int *, std::complex<double> *, int *, int *, int *, int *, std::complex<double> *, int *, int *, int *, int *);
void pdgetri(int *, double *, int *, int *, int *, int *, double*, int *, int *, int *, int *);
void pzgetrf(int *, int *, std::complex<double> *, int *, int *, int *, int *, int *);
void pzgeqpf(int *, int *, std::complex<double> *, int *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *, double *, int *, int*);
void pdgeqpf(int *, int *, double *, int *, int *, int *, int *, double *, double *, int *, int*);
#ifdef __cplusplus
}
#endif


#endif
