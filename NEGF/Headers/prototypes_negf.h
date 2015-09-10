//#include "typedefs.h"
void SetEnergyWeight (std::complex<double> *ene, std::complex<double> *weight, double EF, int *nenergy);
void SetEnergyWeightNoneq (std::complex<double> *ene, std::complex<double> *weight, double EF1, double EF, int *nenergy);
void Run(STATE *, STATE *);

#ifdef __cplusplus
extern "C" {
#endif


void read_LCR();
void read_trans (complex_energy_integral * cei);
void pmo_init();
void lead_bandstructure();
void get_cond_frommatrix_kyz();
void allocate_matrix_soft();

void init_soft (double * vh, double * rho, double * rhocore, double * rhoc, double * rho_tf,
                STATE * states, STATE * states1, double * vnuc, double * vext, double * vxc, double * vh_old,
                double * vxc_old, STATE *states_distribute);

void get_dos(STATE *);
void get_3Ddos(STATE *, double, double, int, int);
void write_matrix_p();

void quench_negf (STATE * states, STATE * states1, STATE *states_distribute, double * vxc, double * vh, double * vnuc, double * vext,
             double * vh_old, double * vxc_old, double * rho, double * rhoc, double * rhocore, double * rho_tf, double * vbias);

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
