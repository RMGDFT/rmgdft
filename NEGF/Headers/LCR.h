/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <complex.h>

#define 	TOLSCF 		1.0E-9


#define 	NUM_PROBE_MAX 	4 


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

NON_LINEAR_THREE_PART lcr[NUM_PROBE_MAX + 1];






void set_energy_weight(complex double *ene, complex double *weight, REAL EF, int *nenergy);
void set_energy_weight_ne(complex double *eneR, complex double *weight, REAL EF1, REAL EF2, int *nenergy);


/* do we need these declarations? I think not... 
 * */
void zgemm(char*, char*, int*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*); 
/* void zgemm(char*, char*, int*, int*, int*, REAL*, REAL*, int*, REAL*, int*, REAL*, REAL*, int*);
 */
void ZCOPY(int*, complex double*, int*, complex double*, int*);
void ZAXPY(int*, complex double*, complex double*, int*, complex double*, int*);
void pzgemm(_fcd, _fcd, int*, int*, int*, complex double*,
        complex double*, int*, int*, int*,
        complex double*, int*, int*, int*,
        complex double*, complex double*, int*, int*, int*);



void Stransfer_f(REAL*, REAL*, REAL*, REAL*, int*);
void mat_part(REAL*, REAL*, REAL*, int);
void part_to_distribute(REAL*, REAL*, int);
void read_LCR();
void read_data_LCR(char *nameL, char *nameC, char *nameR, double *vh,
        double *vxc, double *rho, STATE *states);
void read_data_part(char *name, double *vh, double *vxc, double *rho, int which_part);

void Sigma(REAL *sigma, REAL *HLC, REAL *SLC, REAL eneR, REAL eneI, REAL *Green, int iprobe);

void Sgreen_c_wang(REAL *Htri, REAL *Stri, complex double *sigma_all, int *sigma_idx, REAL eneR, REAL eneI, complex double *Green_C, int nC);

void Sgreen_c(REAL *H00, REAL *S00, complex double *sigma, complex
        double *, complex double ene, complex double *Green_C, int nC);
void distri_fermi(complex double ene, REAL EF, complex double *distri);
void rho_munu(complex double *rho_mn, complex double *green_C, complex double *sigma_L, int nL, int N, int *ni, int iprobe);
void modify_rho(REAL *rho, REAL *rho_old);
void modify_rho_y(REAL *rho, REAL *rho_old);
void read_data_lead(double *vh, double *vxc, double *vh_old, double *vxc_old, double *rho);
void write_data_lead(char *name, double *vh, double *vxc, double *vh_old, double *vxc_old, double *rho);
void get_cond_lead();
void whole_to_tri_real(REAL *A_tri, REAL *Aii, int N, int *ni);

void    get_inverse_block(complex double *Hii, complex double *Gii, int *ipiv, int nn);

void read_data_conductor(char *name, double *vh, double *vxc, double *rho);


//void PZGESV(int*, int*, complex double *, int*, int*, int*, int*,
//                complex double*, int*, int*, int*, int*);

