/****** Compile-time options for real-space code are set here  ******/

/* Number of processors */
#define NPES 27

/* 3D processor grid PE_X*PE_Y*PE_Z must equal NPES */ 
#define PE_X 3
#define PE_Y 3
#define PE_Z 3

/* Gamma point only set to 1, otherwise, 0 */
#define GAMMA_PT 1

/* Number of points in coarse grid in x,y and z directions*/
#define NX_GRID 48
#define NY_GRID 48
#define NZ_GRID 48

/* How many times the fine grid is finer than the coarse grid 
 * All three numbers have to be the same */
#define FG_NX 2
#define FG_NY 2
#define FG_NZ 2

/* Set this to 0 to turn off memory Smart-ALLOCation. (salloc.c, salloc.h)
 * (there is no significant time improvement)*/
#define USE_SALLOC 1

/* Set this to 1 if you want to use finite difference method for calculating
 * derivatives of beta. This is faster since it avoids doing 3 backwards fourier
 * transforms per ion, but it may not be very accurate since the finite diff
 * derivative is done on the coarse grid.
 * Leave this set to 0 unless you know what you are doing */
#define FDIFF_BETA 0

/* This enables pulay mixing. This is experimental, used only for code development */
#define PULAY 0

/* Extra fine timers, may cause some slowdown, but they are generally useful */
#define MD_TIMERS 1
