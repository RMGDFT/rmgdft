/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <complex.h>

#define 	TOLSCF 		1.0E-9


#define 	NUM_PROBE_MAX 	4 
#define 	NUM_SUBSYSTEM_MAX   15

complex double *sigma_all;

#include "cei.h"

struct NON_LINEAR_THREE_PART2
{
 int nenergy_ne;

 complex double *ene_ne;
 complex double *weight_ne;

 rmg_double_t *density_matrix_ne;
 rmg_double_t *density_matrix_ne_tri;

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

 rmg_double_t *density_matrix;

 complex double *ene;
 complex double *weight;


 rmg_double_t bias;
 rmg_double_t EF_old;
 rmg_double_t EF_new;

 int nenergy;

 rmg_double_t *S00;
 rmg_double_t *H00;
 rmg_double_t *S01;
 rmg_double_t *H01;

 rmg_double_t *HCL;
 rmg_double_t *SCL;

 rmg_double_t *Htri;
 rmg_double_t *Stri;
 rmg_double_t *density_matrix_tri;

 rmg_double_t xside;
 rmg_double_t yside;

 rmg_double_t 	x_shift; 
 rmg_double_t 	y_shift; 

 rmg_double_t  *S00_yz;
 rmg_double_t  *H00_yz;
 rmg_double_t  *S01_yz;
 rmg_double_t  *H01_yz;
 rmg_double_t  *HCL_yz;
 rmg_double_t  *SCL_yz;
 rmg_double_t *Htri_yz;
 rmg_double_t *Stri_yz;


NON_LINEAR_THREE_PART2 lcr_ne[NUM_PROBE_MAX];
};
typedef struct  NON_LINEAR_THREE_PART NON_LINEAR_THREE_PART;

NON_LINEAR_THREE_PART lcr[NUM_SUBSYSTEM_MAX];






