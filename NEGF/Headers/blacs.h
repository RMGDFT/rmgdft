/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#if(M_SGI_ORIGIN_MPI)

  #define blacs_pinfo blacs_pinfo_
  #define blacs_setup blacs_setup_
  #define blacs_get blacs_get_
  #define blacs_set blacs_set_
  #define blacs_gridinit blacs_gridinit_
  #define blacs_gridmap blacs_gridmap_
  #define blacs_freebuff blacs_freebuff_
  #define blacs_gridexit blacs_gridexit_
  #define  blacs_abort     blacs_abort_  
  #define  blacs_exit      blacs_exit_
  #define  blacs_gridinfo  blacs_gridinfo_
  #define  blacs_pnum      blacs_pnum_
  #define  blacs_pcoord    blacs_pcoord_
  #define  blacs_barrier   blacs_barrier_
  #define dgerv2d dgerv2d_
  #define dgesd2d dgesd2d_
  #define sgerv2d sgerv2d_
  #define sgesd2d sgesd2d_
  #define igesd2d igesd2d_
  #define igerv2d igerv2d_
  #define igamn2d igamn2d_
#endif

/* initialization */
void Cblacs_pinfo(int*, int*);
void Cblacs_setup(int*, int*);
void Cblacs_get(int, int, int*);
void Cblacs_set(int, int, int*);
void Cblacs_gridinit(int*, char*, int, int);
void Cblacs_gridmap(int*, int*, int, int, int);
/* desctruction */
void Cblacs_freebuff(int,int);
void Cblacs_gridexit(int);
void Cblacs_abort(int,int);
void Cblacs_exit(int);
/* Informational and miscellaneous */
void Cblacs_gridinfo(int, int*, int*, int*, int*);
int Cblacs_pnum(int, int, int);
void Cblacs_pcoord(int,int,int*,int*);
void Cblacs_barrier(int, char*);
void Cdgerv2d(int, int, int, double *, int, int, int);
void Cdgesd2d(int, int, int, double *, int, int, int);
void Cigamn2d(int, char scope[], char top[], int m, int n, int *A, int lda,
              int *ra, int *ca, int rcflag, int rdest, int cdest);


/* initialization */
void blacs_pinfo(int*, int*);
void blacs_setup(int*, int*);
void blacs_get(int*, int*, int*);
void blacs_set(int*, int*, int*);
void blacs_gridinit(int*, char*, int*, int*);
void blacs_gridmap(int*, int*, int*, int*, int*);
/* desctruction */
void blacs_freebuff(int*,int*);
void blacs_gridexit(int*);
void blacs_abort(int*,int*);
void blacs_exit(int*);
/* Informational and miscellaneous */
void blacs_gridinfo(int*, int*, int*, int*, int*);
int blacs_pnum(int*, int*, int*);
void blacs_pcoord(int*,int*,int*,int*);
void blacs_barrier(int, char*);
void dgerv2d(int*, int*, int*, double *, int*, int*, int*);
void dgesd2d(int*, int*, int*, double *, int*, int*, int*);
void sgerv2d(int*, int*, int*, double *, int*, int*, int*);
void sgesd2d(int*, int*, int*, double *, int*, int*, int*);
void igesd2d(int*, int*, int*, int *, int*, int*, int*);
void igerv2d(int*, int*, int*, int*, int*, int*, int*);
void igamn2d(int*, char*, char*, int*, int*, int*, int*,
              int*, int*, int*, int*, int*);


