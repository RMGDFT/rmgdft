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


/* Blacs dimension */
#define DLEN    9

// Maximum number of scalapack groups
#define MAX_SCALAPACK_GROUPS 32
// Maximum number of folded spectrum scalapack instances
#define MAX_FOLDED_SCALAPACKS 32

extern "C" int Csys2blacs_handle(MPI_Comm SysCtxt );
extern "C" MPI_Comm Cblacs2sys_handle (int BlacsCtxt);

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

extern "C" {

int numroc_ (int *, int *, int *, int *, int *);
int INDXG2P (int *, int *, int *, int *, int *);
void descinit_ (int[], int *, int *, int *, int *, int *, int *, int *, int *,
               int *);
void pdgesv_ (int *, int *, double *, int * , int *, int *, int *, double *,
        int *, int *, int *, int *);
void pzgesv_ (int *, int *, double *, int * , int *, int *, int *, double *,
        int *, int *, int *, int *);
void pdgemm_ (char *, char *, int *, int *, int *, double *, double *, int *,
             int *, int *, double *, int *, int *, int *, double *, double *,
             int *, int *, int *);
void pzgemm_ (char *, char *, int *, int *, int *, double *, double *, int *,
             int *, int *, double *, int *, int *, int *, double *, double *,
             int *, int *, int *);
void pdsyev_ (char *, char *, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, int *);
void pdsyevd_ (char *, char *, int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *, int *, int *, int *, int *);
void pzheevd_ (char *, char *, int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *, int *, double *, int *, int *, int *, int *);
void pcheev_ (char *, char *, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, double *, int *, int *);
void pdsygst_(int *, char *, int *, double *, int *, int *, int *, double *, int *,
              int *, int *, double *, int *);
void pdsyngst_(int *, char *, int *, double *, int *, int *, int *, double *, int *,
              int *, int *, double *, double *, int *, int *);

void pzhegst_(int *, char *, int *, double *, int *, int *, int *, double *, int *,
              int *, int *, double *, int *);
void pdpotrf_(char *, int*, double*, int*, int*, int*, int*);
void pzpotrf_(char *, int*, double*, int*, int*, int*, int*);
void pdtrtri_(char *, char *, int*, double*, int*, int*, int*, int*);
void pdsyrk_( char *, char *, int *, int *, double *, double *, int *, int *, int *,
             double *, double *, int *, int *, int *);
void pdlaset_(char *, int *, int *, double *, double *, double *, int *, int *, int *);

void PSPOCON (char *, int *, double *, int *, int *, int *, double *, double *,
              double *, int *, int *, int *, int *);
void PSPOTRF (char *, int *, double *, int *, int *, int *, int *);
void PSPOTRI (char *, int *, double *, int *, int *, int *, int *);
void PSSYGST (int *, char *, int *, double *, int *, int *, int *, double *,
              int *, int *, int *, double *, int *);
void PSTRTRS (char *, char *, char *, int *, int *, double *, int *, int *, int *,
              double *, int *, int *, int *, int *);
void PSSYMM (char *, char *, int *, int *, double *, double *, int *, int *,
             int *, double *, int *, int *, int *, double *, double *, int *,
             int *, int *);

void PSUBDIAG (char *, char *, int, double *, int, double *, int *);
void PDSYGVX(int *, char*, char*, char*, int*, double *, int*, int*, int*, double*, int*, int*,
       int*, double*, double *, int*, int*, double*, int*, int*, double*, double*, double*, int*,
       int*, int*, double*, int*, int*, int*, int*, int*, double*, int*);
void PZHEGVX(int *, char*, char*, char*, int*, double *, int*, int*, int*, double*, int*, int*,
       int*, double*, double *, int*, int*, double*, int*, int*, double*, double*, double*, int*,
       int*, int*, double*, int *, double *, int*, int*, int*, int*, int*, double*, int*);
void PDSYEVX(char*, char*, char*, int*, double *, int*, int*, int*, double*, double*, int*,
       int*, double*, int*, int*, double*, double*, double*, int*,
       int*, int*, double*, int*, int*, int*, int*, int*, double*, int*);
void pdgeadd_(char *, int *, int *, double *, double *, int *, int *, int *, double *,
       double *, int *, int *, int *);               
void pzgeadd_(char *, int *, int *, double *, double *, int *, int *, int *, double *,
       double *, int *, int *, int *);               
void pdtrmm_(char *SIDE, char *UPLO, char *TRANS, char *DIAG, int * M, int *N, double *ALPHA,
             double * A, int *IA, int *JA, int *DESCA, double *B, int *IB, int *JB, int *DESCB);               
void pdtrsm_(char *SIDE, char *UPLO, char *TRANS, char *DIAG, int * M, int *N, double *ALPHA,
             double * A, int *IA, int *JA, int *DESCA, double *B, int *IB, int *JB, int *DESCB);               
void pztrmm_(char *SIDE, char *UPLO, char *TRANS, char *DIAG, int * M, int *N, double *ALPHA,
             double * A, int *IA, int *JA, int *DESCA, double *B, int *IB, int *JB, int *DESCB);               
void pztrsm_(char *SIDE, char *UPLO, char *TRANS, char *DIAG, int * M, int *N, double *ALPHA,
             double * A, int *IA, int *JA, int *DESCA, double *B, int *IB, int *JB, int *DESCB);               


}


#endif
