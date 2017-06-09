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


extern BaseGrid *Rmg_G;
extern TradeImages *Rmg_T;
extern Lattice Rmg_L;
extern MpiQueue *Rmg_Q;

extern Pw *coarse_pwaves, *fine_pwaves, *beta_pwaves;
#include "fft3d.h"

extern struct fft_plan_3d *fft_forward_coarse, *fft_backward_coarse, *fft_forward_fine, *fft_backward_fine;


extern "C"
{
double my_crtc (void);
void thread_barrier_wait(void);
void get_vxc (double * rho, double * rho_oppo, double * rhocore, double * vxc);
void symmetry (int *ibrav, int *s, int *nsym, int *irg, int *irt,
               int *ftau, int *nat, double * tau, int *ityp, int *nks,
               double * xk, double * wk, double * celldm, int *nr1, int *nr2,
               int *nr3, double *a1, double *a2, double *a3, double *omega, int *wflag);
}
extern "C" double get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target, int boundaryflag);
extern "C" void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio);
extern "C" void mix_rho (double * new_rho, double * rho, double *rhocore, int length, int length_x, int length_y, int length_z);
extern "C" void  get_rho_oppo (double * rho, double * rho_oppo);
extern "C" void get_ddd (double *veff);
extern "C" void mix_betaxpsi (int mix);
extern "C" void rmg_lbfgs (void);
extern "C" void write_restart (char *name, double * vh, double *vxc, double *vh_old, 
        double *vxc_old,  double * rho, double *rho_oppo, STATE *states);

extern "C" int init_kpoints (int *mesh, int *is_shift);

template <typename RmgType> void AppCir (RmgType * a, RmgType * b, char * grid);
template <typename RmgType> double AppCil (RmgType * a, RmgType * b, char * grid);
template <typename DataType> double ApplyAOperator (DataType *a, DataType *b, char *grid);
template <typename DataType> double ApplyAOperator (DataType *a, DataType *b, char *grid, BaseGrid *G, TradeImages *T);
template <typename RmgType> void ApplyBOperator (RmgType * a, RmgType * b, char *grid);
template <typename RmgType> void ApplyBOperator (RmgType * a, RmgType * b, char *grid, BaseGrid *G, TradeImages *T);
template <typename DataType> void ApplyGradient (DataType *a, DataType *gx, DataType *gy, DataType *gz, int order, char *grid);
template <typename DataType> void ApplyGradient (DataType *a, DataType *gx, DataType *gy, DataType *gz, int order, char *grid, BaseGrid *G, TradeImages *T);
template <typename DataType> double ApplyLaplacian (DataType *a, DataType *b, int order, char *grid);
template <typename DataType> double ApplyLaplacian (DataType *a, DataType *b, int order, char *grid, BaseGrid *G, TradeImages *T);

void GetVtotPsi (double * vtot_psi, double * vtot, int grid_ratio);
template <typename KpointType>
void MolecularDynamics (Kpoint<KpointType> **Kptr, double * vxc, double * vh, double * vnuc,
             double * rho, double * rho_oppo, double * rhoc, double * rhocore);






template <typename OrbitalType> void GetNewRho(Kpoint<OrbitalType> **Kpts, double *rho);
template <typename OrbitalType> void GetAugRho(Kpoint<OrbitalType> **Kpts, double *rho);
template <typename OrbitalType> void Init (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> void Relax (int steps, double * vxc, double * vh, double * vnuc,
              double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> bool Quench (double * vxc, double * vh, double * vnuc, double * rho,
             double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> bool Scf (double * vxc, double *vxc_correct, double * vh, double *vh_in, double *vh_ext,
          double * vnuc, double * rho, double * rho_oppo, double * rhocore, double * rhoc, int spin_flag,
          int boundaryflag, Kpoint<OrbitalType> **Kptr, std::vector<double>& RMSdV);
template <typename KpointType> void AppNls(Kpoint<KpointType> *kpoint, KpointType *sintR, 
            KpointType *psi, KpointType *nv, KpointType *ns, KpointType *Bns, int first_state, int num_states, bool need_bns);
template <typename KpointType> void AppNls(Kpoint<KpointType> *kpoint, KpointType *sintR, 
            KpointType *psi, KpointType *nv, KpointType *ns, KpointType *Bns, int first_state, int num_states);
template <typename OrbitalType> double EnergyCorrection (Kpoint<OrbitalType> **Kptr,
          double *rho, double *new_rho, double *vh, double *vh_in);



template <typename OrbitalType, typename CalcType> void MgEigState (Kpoint<OrbitalType> *kptr, 
State<OrbitalType> * sp, double * vtot_psi, OrbitalType *nv, OrbitalType *ns, int vcycle);


// Gamma point float version
void CPP_genvpsi (float * psi, float * sg_twovpsi, double * vtot, void * kd,
              double kmag, int dimx, int dimy, int dimz);
// complex float version
void CPP_genvpsi (std::complex<float> * psi, std::complex<float> * sg_twovpsi, double * vtot, void * kd,
              double kmag, int dimx, int dimy, int dimz);
// complex double version
void CPP_genvpsi (std::complex<double> * psi, std::complex<double> * sg_twovpsi, double * vtot, void * kd,
              double kmag, int dimx, int dimy, int dimz);
// Gamma point double version
void CPP_genvpsi (double * psi, double * sg_twovpsi, double * vtot, void * kd,
              double kmag, int dimx, int dimy, int dimz);

void pack_to_complex(double *psi, int nstates, int pbasis);
void pack_to_standard(double *psi, int nstates, int pbasis);
void MixBetaxpsi (int mix, int kpt);
template  <typename OrbitalType> void AppCilrDriver (TradeImages *T, OrbitalType * psi, OrbitalType * a_psi, OrbitalType *b_psi, double *vtot,
    int dimx, int dimy, int dimz, double hx, double hy, double hz, int order);
void MixRho (double * new_rho, double * rho, double *rhocore, double *vh_in, double *vh_out, double *rhoc, std::unordered_map<std::string, InputKey *>& ControlMap, bool reset);

void DiagScalapack(STATE *, int, double *, double*, double *, double *);
void DiagElemental(STATE *, int, double *, double*, double *, double *);
void BandwidthReduction(int num_ions, ION *ions, unsigned int *);
void PermAtoms(int num_ions, ION *ions, unsigned int *);
void GetPermStateIndex(int num_ions, ION *ions, unsigned int *, unsigned int *, unsigned int *);
void ReadPermInfo(char *, unsigned int *);
void WritePermInfo(char *, unsigned int *);
void InitPe4image();
void InitPe4kpspin();


extern "C" void app_cilr_driver (double * psi, double * a_psi, double *b_psi, double *vtot_eig_s,
    int dimx, int dimy, int dimz, double hx, double hy, double hz, int order);

template  <typename OrbitalType> double AppCilrFourth (OrbitalType *psi, OrbitalType *a_psi, OrbitalType *b_psi, double *vtot, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);
template  <typename OrbitalType> double AppCilrSixth (OrbitalType *psi, OrbitalType *a_psi, OrbitalType *b_psi, double *vtot, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);
template <typename OrbitalType, typename CalcType>
void PotentialAcceleration(Kpoint<OrbitalType> *kptr, State<OrbitalType> *sp, double *vtot_psi, double *nvtot_psi, CalcType *tmp_psi_t, OrbitalType *saved_psi);
// Print function
void RmgPrintTimings(MPI_Comm comm, const char *outfile, int steps, int num_ions_loc, int override_rank);
void RmgPrintTimings(MPI_Comm comm, const char *outfile, int steps, int num_ions_loc);
template <typename KpointType>
void ReinitIonicPotentials (Kpoint<KpointType> **kptr, double * vnuc, double * rhocore, double * rhoc);
template <typename KpointType>
void GetNlop (Kpoint<KpointType> **Kptr);
template <typename KpointType>
void GetWeight (Kpoint<KpointType> **Kptr);
template <typename KpointType>
void GetWeightLocal (Kpoint<KpointType> **Kptr);

template <typename KpointType>
void AssignWeight (Kpoint<KpointType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, double * rtptr, KpointType *Bweight, KpointType *Nlweight);
template <typename KpointType>
void Betaxpsi (Kpoint<KpointType> *kptr, int, int, KpointType *, KpointType *);
template <typename RmgType>
void AppCirDriverBeta (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order);
template <typename DataType>
double ApplyAOperator (Lattice *L, TradeImages *T, DataType *a, DataType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);
template <typename RmgType>
void ApplyBOperator (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order);
void LoadUpf(SPECIES *sp);
extern "C" void LoadUpf_C(SPECIES *sp);
extern "C" bool verify( char *tagname, const void *optvalue );
void ReadPseudo(int nspecies, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
template <typename KpointType>
void OutputEigenvalues (Kpoint<KpointType> **Kptr, int ikbs, int iscf);
template <typename KpointType>
void ReadData (char *name, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
void GetOppositeEigvals (Kpoint<KpointType> **Kptr);
template <typename KpointType>
void GetOppositeOccupancies (Kpoint<KpointType> **Kptr);
template <typename StateType>
void LcaoGetPsi (State<StateType> * states);
template <typename StateType>
void LcaoGetAwave (StateType *psi, ION *iptr, int awave_idx, int l, int m, double coeff);
void LcaoGetRho (double * arho_f);
template <typename KpointType>
void GetTe (double * rho, double * rho_oppo, double * rhocore, double * rhoc, double * vh, double * vxc, Kpoint<KpointType> ** Kptr , int ii_flag);
template <typename KpointType>
void WriteRestart (char *name, double * vh, double * rho, double * rho_oppo, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
void WriteBGW_Wfng (int kpt, Kpoint<KpointType> * Kptr);
void WriteBGW_Rhog (double * rho, double * rho_oppo);
template <typename KpointType>
void WriteBGW_VxcEig (int kpt, double * vxc, Kpoint<KpointType> * Kptr);
template <typename KpointType>
void WriteData (int fhand, double * vh, double * rho, double * rho_oppo, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
double Fill (Kpoint<KpointType> **Kptr, double width, double nel, double mix, int num_st, int occ_flag, int mp_order);
template <typename KpointType>
void OutputBandPlot(Kpoint<KpointType> ** Kptr);
int GetRecommendedThreadNumber(int nthreads, int npes, int thispe, MPI_Comm comm);
void InitHybridModel(int nthreads, int npes, int thispe, MPI_Comm comm);
void ReadCommon(int argc, char *argv[], char *cfile, CONTROL& cont, PE_CONTROL& pecont, std::unordered_map<std::string, InputKey *>& Map);
void AutoSet(CONTROL& cont, PE_CONTROL& pecont, std::unordered_map<std::string, InputKey *>& Map);
void ReadInit(char *meta, CONTROL& lc, PE_CONTROL& pelc, std::unordered_map<std::string, InputKey *>& InputMap);
void InitIo (int argc, char **argv, std::unordered_map<std::string, InputKey *>& Map);
bool Verify(const std::string& KeyName, const std::string& KeyVal, const std::unordered_map<std::string, InputKey *>& Map);
bool Verify(const std::string& KeyName, const char *keyval, const std::unordered_map<std::string, InputKey *>& Map);
bool Verify(const std::string& KeyName, const bool& KeyVal, const std::unordered_map<std::string, InputKey *>& Map);
void ReadDynamics(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadRmgAtoms(char *cfile, std::set<std::string>& SpeciesTypes, std::list<std::string>& Species, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadForces(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadVelocities(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadTFAtoms(char *cfile, std::set<std::string>& SpeciesTypes, std::list<std::string>& IonSpecies, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);



template <typename OrbitalType> 
    void Force (double * rho, double * rho_oppo, double * rhoc, double * vh, double *vh_in,
        double * vxc, double *vxc_in, double * vnuc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> 
    void Nlforce (double * , Kpoint<OrbitalType> **Kptr, double *force);


template <typename OrbitalType> void PartialBetaxpsi (int ion, fftw_plan p2, double * newsintR_x, double * newsintR_y,
                       double * newsintR_z, double * newsintI_x, double * newsintI_y, double * newsintI_z,
                       ION * iptr,Kpoint<OrbitalType> **Kptr);

template <typename OrbitalType> void GetGamma (double * gammaR, int ion, int nh , Kpoint<OrbitalType> **Kptr);


template <typename OrbitalType> void PartialGamma (
                    int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
                    Kpoint<OrbitalType> **kptr, int state_start, int state_end, 
                    OrbitalType *sint_derx, OrbitalType *sint_dery, OrbitalType *sint_derz);



template <typename OrbitalType> void AssignDerweight (Kpoint<OrbitalType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, OrbitalType
*rtptr);
        

void ReadKpoints(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
int ReadKpointsBandstructure(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
void ReadOrbitals(char *cfile, STATE  *states, ION *ions,  MPI_Comm comm, unsigned int *);
void ReadBranchON(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap);
template <typename KpointType>
void BandStructure(Kpoint<KpointType> ** Kptr, double *vh, double *vxc, double *vnuc);
void GetPrimeFactors(std::vector<int>& factors, int val, int stop);
void SetupGrids(int npes, int& NX_GRID, int& NY_GRID, int& NZ_GRID, double *celldm, double h, PE_CONTROL& pelc);
void SetupProcessorGrid(int npes, int NX_GRID, int NY_GRID, int NZ_GRID, PE_CONTROL& pelc);
void SetupWavefunctionGrid(int npes, int& NX_GRID, int& NY_GRID, int& NZ_GRID, double *celldm, double h);
extern "C" int FilenameIncrement(char *pathname);

template <typename RmgType>
void CPP_app_smooth (RmgType * f, RmgType * work, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_app_smooth1 (RmgType * f, RmgType * work, int dimx, int dimy, int dimz);

int CountAtomicOrbitals(void);
void InitPseudo (std::unordered_map<std::string, InputKey *>& ControlMap);
void InitQfunct (std::unordered_map<std::string, InputKey *>& ControlMap);
void FindPhase (int nlxdim, int nlydim, int nlzdim, double * nlcdrs, std::complex<double> *phase);
void VhPfft(double *rho, double *rhoc, double *vh);
double VhDriver(double *rho, double *rhoc, double *vh, double *vh_ext, double rms_target);
void BroydenPotential(double *rho, double *new_rho, double *rhoc, double *vh_in, double *vh_out, int max_iter, bool reset);
void output_force(double *force, char *desc);
int Radius2grid (double radius, double mingrid_spacing);
void Lforce (double * rho, double * vh, double *force);
void Nlccforce (double * rho, double * vxc, double *force_nlcc);
void CorrectForces (double *vh, double *vh_in, double *vxc, double *vxc_in, double *force);

void FindFftwPhaseLocalpp (int nlxdim, int nlydim, int nlzdim,double * nlcdrs, std::complex<double> *phase_fft, int level);
void InitLocalForward();
void InitLocalBackward(double *, double *, double *);
void GetPhase (ION *, std::complex<double> *);
void PrintSums(double *, int, char*);
void InitWeightOne (SPECIES * sp, fftw_complex * rtptr, std::complex<double> *phaseptr, int ip, int l, int m, 
        fftw_plan p1, fftw_plan p2);
void InitWeight(void);
void InitDelocalizedWeight (void);
double CubicHarmonic(int L, int M, double *r);
double Ylm(int L, int M, double *r);
void InitClebschGordan (int lmax, double *ap, int *lpx, int *lpl);
void GetQI(void);
void GetPhaseSpecies(SPECIES *sp, std::complex<double> *phaseptr);
void PackGftoc(int, int, int, int, int, int, std::complex<double> *, std::complex<double> *);
double IonIonEnergy_Ewald();
void IIforce(double *);
void InitDelocalizedLocalpp(double *vlocpp_r);
template <typename DataType> void AppGradPfft (DataType *a, DataType *gx, DataType *gy, DataType *gz, char *grid);

#endif
#endif

#if !(defined(_WIN32) || defined(_WIN64))
    #define rmg_printf( message... ) \
         fprintf( ct.logfile, message )
#else
    #define rmg_printf( message, ... ) \
         fprintf( ct.logfile, message, __VA_ARGS__ )
#endif


