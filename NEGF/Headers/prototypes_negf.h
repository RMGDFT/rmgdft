//#include "typedefs.h"
void SetEnergyWeight (std::complex<double> *ene, std::complex<double> *weight, double EF, int *nenergy);
void SetEnergyWeightNoneq (std::complex<double> *ene, std::complex<double> *weight, double EF1, double EF, int *nenergy);
void Run(STATE *, STATE *);
void QuenchNegf (STATE * states, STATE * states1, STATE *states_distribute, double * vxc, double * vh, double * vnuc, double * vext,
             double * vh_old, double * vxc_old, double * rho, double * rhoc, double * rhocore, double * rho_tf, double * vbias);
void InitNegf (double * vh, double * rho, double * rhocore, double * rhoc, double * rho_tf,
                STATE * states, STATE * states1, double * vnuc, double * vext, double * vxc, double * vh_old,
                double * vxc_old, STATE *states_distribute);

#ifdef __cplusplus
extern "C" {
#endif

void read_orbital(STATE*);
void interpolation_orbit(STATE*);
void init_state_distribute(STATE *, STATE*);
void scale_orbital(STATE *, STATE*);
void allocate_matrix_LCR();
void read_lead_matrix();
void read_potrho_LCR (double *, double *, double *);
void read_rho_and_pot (char*, double *, double *, double *, double *, double *);
void plane_average_rho(double *); 
void write_rho_x (double *, char*);  // information about vh_init is recorded in zvec array
void init_comp (double *);
void init_psp_soft ();
void init_ext(double *, double, double, double, double);

void zero_lead_image(double*);
void setback_corner_matrix_H();
void setback_corner_matrix_S();
void row_to_tri_p(double *, double *, int, int*);
void sigma_all_energy_point(DoubleC *, double, double);

void scf (DoubleC * sigma_all, STATE * states, STATE * states_distribute, double *vxc,
        double *vh, double *vnuc, double *vext, double *rho, double *rhoc, double *rhocore, double *rho_tf,
        double * vxc_old, double * vh_old, double * vbias, int *CONVERGENCE);


void read_LCR();
void read_trans (complex_energy_integral * cei);
void pmo_init();
void lead_bandstructure();
void get_cond_frommatrix_kyz();
void allocate_matrix_soft();


void get_dos(STATE *);
void get_3Ddos(STATE *, double, double, int, int);
void write_matrix_p();


void apply_potential_drop(double *);


void gauleg(double , double, double *, double *, int);
void writeout_matrix_p();
void write_data_lead(char *, double *, double *, double *, double *, double *);
void write_data_negf(char *, double *,double *,double *,double *,double *,double *, STATE *);
void set_energy_weight (DoubleC *ene, DoubleC *weight, double EF, int *nenergy);
void set_energy_weight_ne (DoubleC *ene, DoubleC *weight, double EF1, double EF, int *nenergy);

#ifdef __cplusplus
}
#endif


