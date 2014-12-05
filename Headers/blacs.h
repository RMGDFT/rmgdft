#ifndef RMG_blacs_h
#define RMG_blacs_h 1

#if __cplusplus
extern "C" {
#endif
 
/* initialization */
void Cblacs_pinfo (int *, int *);
void Cblacs_setup (int *, int *);
void Cblacs_get (int, int, int *);
void Cblacs_set (int, int, int *);
void Cblacs_gridinit (int *, char *, int, int);
void Cblacs_gridmap (int *, int *, int, int, int);
void Cblacs_gridinfo (int, int *, int *, int *, int *);
/* desctruction */
void Cblacs_freebuff (int, int);
void Cblacs_gridexit (int);
void Cblacs_abort (int, int);
void Cblacs_exit (int);
/* Informational and miscellaneous */
void Cblacs_gridinfo (int, int *, int *, int *, int *);
int Cblacs_pnum (int, int, int);
void Cblacs_pcoord (int, int, int *, int *);
void Cblacs_barrier (int, char *);
void Cdgerv2d (int, int, int, double *, int, int, int);
void Cdgesd2d (int, int, int, double *, int, int, int);
void Cigamn2d (int, char scope[], char top[], int m, int n, int *A, int lda,
               int *ra, int *ca, int rcflag, int rdest, int cdest);



/* initialization */
void blacs_pinfo (int *, int *);
void blacs_setup (int *, int *);
void blacs_get (int *, int *, int *);
void blacs_set (int *, int *, int *);
void blacs_gridinit (int *, char *, int *, int *);
void blacs_gridmap (int *, int *, int *, int *, int *);
/* desctruction */
void blacs_freebuff (int *, int *);
void blacs_gridexit (int *);
void blacs_abort (int *, int *);
void blacs_exit (int *);
/* Informational and miscellaneous */
void blacs_gridinfo (int *, int *, int *, int *, int *);
int blacs_pnum (int *, int *, int *);
void blacs_pcoord (int *, int *, int *, int *);
void blacs_barrier (int, char *);
void dgerv2d (int *, int *, int *, double *, int *, int *, int *);
void dgesd2d (int *, int *, int *, double *, int *, int *, int *);
void sgerv2d (int *, int *, int *, double *, int *, int *, int *);
void sgesd2d (int *, int *, int *, double *, int *, int *, int *);
void igesd2d (int *, int *, int *, int *, int *, int *, int *);
void igerv2d (int *, int *, int *, int *, int *, int *, int *);
void igamn2d (int *, char *, char *, int *, int *, int *, int *,
              int *, int *, int *, int *, int *);

#if __cplusplus
}
#endif

#endif
