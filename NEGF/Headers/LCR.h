/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <complex.h>

#define 	TOLSCF 		1.0E-9


#define 	NUM_PROBE_MAX 	4 
#define 	NUM_SUBSYSTEM_MAX   15

extern double _Complex *sigma_all;

struct NON_LINEAR_THREE_PART2
{
 int nenergy_ne;

 double _Complex *ene_ne;
 double _Complex *weight_ne;

 double *density_matrix_ne;
 double *density_matrix_ne_tri;

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

 int num_states;
 int num_ions;
 int state_begin;
 int state_middle;
 int state_end;
 int ion_begin;
 char name[MAX_PATH];
 char lead_name[MAX_PATH];

 double *density_matrix;

 double _Complex *ene;
 double _Complex *weight;


 double bias;
 double EF_old;
 double EF_new;

 int nenergy;

 double *S00;
 double *H00;
 double *S01;
 double *H01;

 double *HCL;
 double *SCL;

 double *Htri;
 double *Stri;
 double *density_matrix_tri;

 double xside;
 double yside;

 double 	x_shift; 
 double 	y_shift; 

 double  *S00_yz;
 double  *H00_yz;
 double  *S01_yz;
 double  *H01_yz;
 double  *HCL_yz;
 double  *SCL_yz;
 double *Htri_yz;
 double *Stri_yz;


NON_LINEAR_THREE_PART2 lcr_ne[NUM_PROBE_MAX];
};
typedef struct  NON_LINEAR_THREE_PART NON_LINEAR_THREE_PART;

extern NON_LINEAR_THREE_PART lcr[NUM_SUBSYSTEM_MAX];



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
    int probe_noneq;
    int energy_point_insert;
};
typedef struct  complex_energy_integral complex_energy_integral;

extern complex_energy_integral cei;
