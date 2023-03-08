#ifndef RMG_transition_h
#define RMG_transition_h

#if __cplusplus
#include <complex>
#include <set>
#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "Pw.h"
#include "typedefs.h"
#include "LaplacianCoeff.h"
#include "PulayMixing.h"
#include "Symmetry.h"
#include "rmgfiles.h"


extern PulayMixing *Pulay_rho;
extern PulayMixing *Pulay_ldau;
extern PulayMixing *Pulay_orbital;
extern BaseGrid *Rmg_G;
extern BaseGrid *Rmg_halfgrid;
extern TradeImages *Rmg_T;
extern Lattice Rmg_L;
extern MpiQueue *Rmg_Q;
extern Symmetry *Rmg_Symm;

extern Pw *coarse_pwaves, *fine_pwaves, *beta_pwaves, *ewald_pwaves, *half_pwaves;


void CheckSetDefault();
void DistributeToGlobal(double *vtot_c, double *vtot_global);
void OptimizeFdCoeff();
void check_tests();
void write_rho_z(double *, char *);
void ProgressTag(double, double);
void ExxProgressTag(double, double);
double my_crtc (void);
void get_vxc (double * rho, double * rho_oppo, double * rhocore, double * vxc);
void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio);
void mix_rho (double * new_rho, double * rho, double *rhocore, int length, int length_x, int length_y, int length_z);
void  get_rho_oppo (double * rho, double * rho_oppo);
void get_ddd (double *veff, double *vxc, bool ddd0_flag);

void lbfgs_init(int num_coeff);
void write_restart (char *name, double * vh, double *vxc, double *vh_old, 
        double *vxc_old,  double * rho, double *rho_oppo, STATE *states);

int init_kpoints (int *mesh, int *is_shift);

template <typename OrbitalType> void STM_calc (Kpoint<OrbitalType> **Kptr, double *rho, std::vector<double> bias_list, std::vector<double>
height_list);
template <typename DataType> void OutputCubeFile(DataType *a, int grid, std::string filename);
template <typename DataType> double ApplyAOperator (DataType *a, DataType *b);
template <typename DataType> double ApplyAOperator (DataType *a, DataType *b, double *kvec);
template <typename DataType> double ApplyAOperator (DataType *a, DataType *b, int, int, int, double, double, double, int, double *kvec);
template <typename DataType> void ApplyGradient (DataType *a, DataType *gx, DataType *gy, DataType *gz, int order, const char *grid);
template <typename DataType> void SumGradientKvec (DataType *a, DataType *b, double *kvec, const char *grid);
template <typename DataType> void ApplyGradient (DataType *a, DataType *gx, DataType *gy, DataType *gz, int order, const char *grid, BaseGrid *G, TradeImages *T);
template <typename DataType> void ApplyGradient (DataType *a, DataType *gx, DataType *gy, DataType *gz, int dimx, int dimy, int dimz, int order);
template <typename DataType> double ApplyLaplacian (DataType *a, DataType *b, int order, const char *grid);
template <typename DataType> double ApplyLaplacian (DataType *a, DataType *b, int order, const char *grid, BaseGrid *G, TradeImages *T);

void GetVtotPsi (double * vtot_psi, double * vtot, int grid_ratio);


// Gamma point float version
void CPP_genvpsi (float * psi, float * sg_twovpsi, double * vtot, double kmag, int dimx, int dimy, int dimz);
// complex float version
void CPP_genvpsi (std::complex<float> * psi, std::complex<float> * sg_twovpsi, double * vtot, double kmag, int dimx, int dimy, int dimz);
// complex double version
void CPP_genvpsi (std::complex<double> * psi, std::complex<double> * sg_twovpsi, double * vtot, double kmag, int dimx, int dimy, int dimz);
// Gamma point double version
void CPP_genvpsi (double * psi, double * sg_twovpsi, double * vtot, double kmag, int dimx, int dimy, int dimz);

void MixRho (double * new_rho, double * rho, double *rhocore, double *vh_in, double *vh_out, double *rhoc, std::unordered_map<std::string, InputKey *>& ControlMap, bool reset);

void MixLdaU (int ns_size, double * new_ns_occ, double * ns_occ, std::unordered_map<std::string, InputKey *>& ControlMap, bool reset);

void DiagScalapack(STATE *, int, double *, double*);
void DiagGpu(STATE *, int, double *, double*, double *, double *, double *);
void DiagElemental(STATE *, int, double *, double*, double *, double *);
void BandwidthReduction(int num_ions, std::vector<ION> &ions, unsigned int *);
void PermAtoms(int num_ions, std::vector<ION> &ions, unsigned int *);
void GetPermStateIndex(int num_ions, std::vector<ION> &ions, unsigned int *, unsigned int *, unsigned int *);
void ReadPermInfo(char *, unsigned int *);
void WritePermInfo(char *, unsigned int *);
void InitPe4image();
void InitPe4kpspin();

//void RmgDft3(std::vector<ION> &ion);
void RmgDftd3(double *disp, double *grads, double *stress, int version);


template  <typename OrbitalType> double AppCilrFourth (OrbitalType *psi, OrbitalType *a_psi, OrbitalType *b_psi, double *vtot, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);
template  <typename OrbitalType> double AppCilrSixth (OrbitalType *psi, OrbitalType *a_psi, OrbitalType *b_psi, double *vtot, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

// Print function
void RmgPrintTimings(MPI_Comm comm, const char *outfile, int steps, int num_ions_loc, int override_rank);
void RmgPrintTimings(MPI_Comm comm, const char *outfile, int steps, int num_ions_loc);
template <typename KpointType>
void ReinitIonicPotentials (Kpoint<KpointType> **kptr, double * vnuc, double * rhocore, double * rhoc);

template <typename KpointType>
void AssignWeight (Kpoint<KpointType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, KpointType *Nlweight);
template <typename KpointType>
void Betaxpsi (Kpoint<KpointType> *kptr, int, int, KpointType *);
template <typename KpointType>
void LdaplusUxpsi (Kpoint<KpointType> *kptr, int, int, KpointType *);
template <typename DataType>
double ApplyAOperator (Lattice *L, TradeImages *T, DataType *a, DataType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);
void LoadPseudo(SPECIES *sp);
void LoadUpfPseudo(SPECIES *sp);
void LoadXmlPseudo(SPECIES *sp);
double * UPF_str_to_double_array(std::string str, int max_count, int start);
extern "C" void LoadUpf_C(SPECIES *sp);
extern "C" bool verify( char *tagname, const void *optvalue );
void ReadPseudo(int nspecies, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);

int GetRecommendedThreadNumber(int nthreads, int npes, int thispe, MPI_Comm comm);
void InitHybridModel(int omp_nthreads, int mg_threads, int npes, int thispe, MPI_Comm comm);
void ReadCommon(char *cfile, CONTROL& cont, PE_CONTROL& pecont, std::unordered_map<std::string, InputKey *>& Map);
void WriteInputOptions(std::unordered_map<std::string, InputKey *>& Map, std::string type);
void AutoSet(CONTROL& cont, PE_CONTROL& pecont, std::unordered_map<std::string, InputKey *>& Map);
void ReadInit(char *meta, CONTROL& lc, PE_CONTROL& pelc, std::unordered_map<std::string, InputKey *>& InputMap);
void InitIo (int argc, char **argv, std::unordered_map<std::string, InputKey *>& Map);
bool Verify(const std::string& KeyName, const std::string& KeyVal, const std::unordered_map<std::string, InputKey *>& Map);
bool Verify(const std::string& KeyName, const char *keyval, const std::unordered_map<std::string, InputKey *>& Map);
bool Verify(const std::string& KeyName, const bool& KeyVal, const std::unordered_map<std::string, InputKey *>& Map);
extern "C" bool verify_boolean(char *tagname, const void *optvalue );
void ReadDynamics(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadRmgAtoms(char *cfile, std::set<std::string>& SpeciesTypes, std::list<std::string>& Species, std::vector<ION> &Atoms,
        CONTROL& lc);
void ReadForces(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadVelocities(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadTFAtoms(char *cfile, std::set<std::string>& SpeciesTypes, std::list<std::string>& IonSpecies, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void CheckShutdown(void);



void ReadKpoints(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
int ReadKpointsBandstructure(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadOrbitals(char *cfile, STATE  *states, std::vector<ION> &ions,  MPI_Comm comm, unsigned int *);
void ReadBranchON(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void GetPrimeFactors(std::vector<int>& factors, int val, int stop);
void SetupGrids(int npes, int& NX_GRID, int& NY_GRID, int& NZ_GRID, double *celldm, double h, PE_CONTROL& pelc);
void SetupProcessorGrid(int npes, int NX_GRID, int NY_GRID, int NZ_GRID, PE_CONTROL& pelc);
void SetupWavefunctionGrid(int& NX_GRID, int& NY_GRID, int& NZ_GRID, double *celldm, double h);
extern "C" int FilenameIncrement(char *pathname);

int CountAtomicOrbitals(void);
void InitSpinOrbit();
void FindPhase (SPECIES *sp, int nlxdim, int nlydim, int nlzdim, double * nlcdrs, std::complex<double> *phase);
void FindPhaseKpoint (double *kvec, int nlxdim, int nlydim, int nlzdim, double * nlcdrs, std::complex<double>* phase_fftw, bool localize);
void FindNlcrds (double *xtal, double *nlcrds,
                 int nlxdim, int nlydim, int nlzdim,
                 int nlxcstart, int nlycstart, int nlzcstart);
void VhPfft(double *rho, double *rhoc, double *vh);
double VhDriver(double *rho, double *rhoc, double *vh, double *vh_ext, double rms_target);
void BroydenPotential(double *rho, double *new_rho, double *rhoc, double *vh_in, double *vh_out, int max_iter, bool reset);
void output_force(double *force, char *desc);
int Radius2grid (const double radius, const double mingrid_spacing, const int ibrav, const bool is_localized);
void Lforce (double * rho, double * vh, double *force);
void Nlccforce (double * rho, double * vxc, double *force_nlcc);
void CorrectForces (double *vh, double *vh_in, double *vxc, double *vxc_in, double *force);

void FindFftwPhaseLocalpp (int nlxdim, int nlydim, int nlzdim,double * nlcdrs, std::complex<double> *phase_fft, int level);
void GetPhase (ION *, std::complex<double> *);
void PrintSums(double *, int, char*);
void InitWeightOne (SPECIES * sp, fftw_complex * rtptr, std::complex<double> *phaseptr, int ip, int l, int m);
double CubicHarmonic(int L, int M, double *r);
double Ylm(int L, int M, double *r);
void InitClebschGordan (int lmax, double *ap, int *lpx, int *lpl);
void GetQI(void);
void GetPhaseSpecies(SPECIES *sp, std::complex<double> *phaseptr);
void PackGftoc(int, int, int, int, int, int, std::complex<double> *, std::complex<double> *);
double IonIonEnergy_Ewald();
void IIforce(double *);
void InitDelocalizedLocalpp(double *vlocpp_r);
template <typename DataType> void AppGradPfft (DataType *a, DataType *gx, DataType *gy, DataType *gz, const char *grid);

void SetLaplacian();
void WriteHeader (void);
template <typename T> void AppExx(Kpoint<T> *kptr, T *psi, int N, T *vexx, T* nv);
void DeviceSynchronize(void);
void Precond_drho(double *);
template <typename T> void Write_Wfs_forWannier(int kpt_global, Kpoint<T> *kptr, std::vector<bool> exclude_bands, std::string wavefule);
double GetPlanarAnisotropy(double *density);
void GetFdFactor(int kidx);

#endif
#endif

#if !(defined(_WIN32) || defined(_WIN64))
    #define rmg_printf( message... ) \
         fprintf( ct.logfile, message )
#else
    #define rmg_printf( message, ... ) \
         fprintf( ct.logfile, message, __VA_ARGS__ )
#endif


