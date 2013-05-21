/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <complex.h>

#define 	TOLSCF 		1.0E-9


#define 	NUM_PROBE_MAX 	4 
#define 	NUM_SUBSYSTEM_MAX   15


struct complex_energy_integral
{
    double EB;
    double KT;
    double GAMMA;
    double DELTA2;
    double DELTA;
    int ncircle;
    int nmax_gq1;
    int nmax_gq2;
    int num_probe;
	int *probe_in_block;
    int num_subsystem;
	int *subsystem_idx;
    int num_probe_window;
    int *probe_window_start;
    int *probe_window_end;
    int num_dos_window;
    int *dos_window_start;
    int *dos_window_end;

    int Npulaysave;
    int Npulayrefresh;
    double pulaymix;
};
typedef struct  complex_energy_integral complex_energy_integral;

extern complex_energy_integral cei;


struct NON_LINEAR_THREE_PART2
{
 int nenergy_ne;

 complex double *ene_ne;
 complex double *weight_ne;

 REAL *density_matrix_ne;
 REAL *density_matrix_ne_tri;

};
typedef struct  NON_LINEAR_THREE_PART2 NON_LINEAR_THREE_PART2;


struct NON_LINEAR_THREE_PART
{
 short int x0;
 short int x1;
 short int x2;
 short int y0;
 short int y1;
 short int y2;
 short int z0;
 short int z1;
 short int z2;

 short int NX_GRID;
 short int NY_GRID;
 short int NZ_GRID;

 short int num_states;
 short int num_ions;
 short int state_begin;
 short int state_middle;
 short int state_end;
 short int ion_begin;
 char name[MAX_PATH];
 char lead_name[MAX_PATH];

 REAL *density_matrix;

 complex double *ene;
 complex double *weight;


 REAL bias;
 REAL EF_old;
 REAL EF_new;

 int nenergy;

 REAL *S00;
 REAL *H00;
 REAL *S01;
 REAL *H01;

 REAL *HCL;
 REAL *SCL;

 REAL *Htri;
 REAL *Stri;
 REAL *density_matrix_tri;

 REAL xside;
 REAL yside;

 REAL 	x_shift; 
 REAL 	y_shift; 

 REAL  *S00_yz;
 REAL  *H00_yz;
 REAL  *S01_yz;
 REAL  *H01_yz;
 REAL  *HCL_yz;
 REAL  *SCL_yz;
 REAL *Htri_yz;
 REAL *Stri_yz;


NON_LINEAR_THREE_PART2 lcr_ne[NUM_PROBE_MAX];
};
typedef struct  NON_LINEAR_THREE_PART NON_LINEAR_THREE_PART;

NON_LINEAR_THREE_PART lcr[NUM_SUBSYSTEM_MAX];






