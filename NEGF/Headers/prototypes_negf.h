//#include "typedefs.h"
void SetEnergyWeight (std::complex<double> *ene, std::complex<double> *weight, double EF, int *nenergy);
void SetEnergyWeightNoneq (std::complex<double> *ene, std::complex<double> *weight, double EF1, double EF, int *nenergy);
void QuenchNegf (STATE * states, STATE * states1,  double * vxc, double * vh, double * vnuc, double * vext,
             double * vh_old, double * vxc_old, double * rho, double * rhoc, double * rhocore, double * rho_tf, double * vbias);

void KrylovSigma(int n, std::complex<double> *H00, std::complex<double> *H10, 
        std::complex<double> *H01, std::complex<double> *sigma, double lamda_min);
#ifdef __cplusplus

#include <unordered_map>
#include "InputKey.h"

void InitNegf (double * vh, double * rho, double * rhocore, double * rhoc, double * rho_tf,
                STATE * states, STATE * states1, double * vnuc, double * vext, double * vxc, double * vh_old,
                double * vxc_old, std::unordered_map<std::string, InputKey *>& ControlMap);
void Run(STATE *, STATE *, std::unordered_map<std::string, InputKey *>& ControlMap);


extern "C" {
#endif

void read_orbital(STATE*);
void interpolation_orbit(STATE*);
void scale_orbital(STATE *);
void allocate_matrix_LCR();
void read_lead_matrix();
void read_potrho_LCR (double *, double *, double *);
void read_rho_and_pot (char*, double *, double *, double *, double *, double *);
void plane_average_rho(double *); 
void write_rho_x (double *, char*);  // information about vh_init is recorded in zvec array
void init_comp (double *);
void init_ext(double *, double, double, double, double);

void zero_lead_image(double*);
void setback_corner_matrix_H();
void setback_corner_matrix_S();
void row_to_tri_p(double *, double *, int, int*);
void sigma_all_energy_point(DoubleC *, double, double);

void ScfNegf (DoubleC * sigma_all, STATE * states, double *vxc,
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

void get_ddd_update (double *);
void HijUpdate (STATE *, double *, double*);
void find_fermi (DoubleC *);
void charge_density_matrix_p (DoubleC*);
void get_new_rho_local (STATE*, double *);
void get_new_rho_soft (STATE*, double *);
void modify_rho (double *, double *);
void get_vh_negf (double*, double*, double*,int, int, int, double);
void tri_to_row (double *, double *, int, int *);


#ifdef __cplusplus
}
#endif


