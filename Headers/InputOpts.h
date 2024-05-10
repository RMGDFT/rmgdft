#ifndef RMG_InputOpts_H
#define RMG_InputOpts_H 1

#include <unordered_map>
#include <const.h>

static std::unordered_map<std::string, int> bravais_lattice_type = {
        {"None", 0},
        {"Cubic Primitive", 1},
        {"Cubic Face Centered", 2},
        {"Cubic Body Centered", 3},
        {"Hexagonal Primitive", 4},
//        {"Hexagonal Rhombohedral (Trigonal)", 5},
        {"Tetragonal Primitive", 6},
//        {"Tetragonal Body Centered", 7},
        {"Monoclinic Primitive", 12},
        {"Orthorhombic Primitive", 8},
//        {"Orthorhombic Base Centered", 9},
//        {"Orthorhombic Body Centered", 10},
//        {"Orthorhombic Face Centered", 11},
//        {"Monoclinic Base Centered", 13},
        {"Triclinic Primitive", 14}};

static std::unordered_map<std::string, int> tddft_mode = {
        {"electric field", EFIELD},
        {"point charge", POINT_CHARGE}};

static std::unordered_map<std::string, int> atomic_orbital_type = {
        {"localized", LOCALIZED},
        {"delocalized", DELOCALIZED}};

static std::unordered_map<std::string, int> internal_pseudo_type = {
        {"ultrasoft", ULTRASOFT_GBRV},
        {"sg15", NORM_CONSERVING_SG15},
        {"nc_accuracy", NORM_CONSERVING_ACCURACY},
        {"nc_standard", NORM_CONSERVING_STANDARD},
        {"all_electron", ALL_ELECTRON}};

static std::unordered_map<std::string, int> energy_output_units = {
        {"Hartrees", 0},
        {"Rydbergs", 1}};

static std::unordered_map<std::string, int> drho_precond_type = {
        {"Resta", PRECOND_RESTA},
        {"Kerker", PRECOND_KERKER}};

static std::unordered_map<std::string, int> crds_units = {
        {"Bohr", 0},
        {"Angstrom", 1}};

static std::unordered_map<std::string, int> vdw_corr = {
        {"None", 0},
        {"DFT-D2", 1},
        {"Grimme-D2", 1},
        {"DFT-D3", 2},
        };

static std::unordered_map<std::string, int> lattice_units = {
        {"Bohr", 0},
        {"Angstrom", 1},
        {"Alat", 2}};

static std::unordered_map<std::string, int> charge_mixing_type = {
        {"Linear", 0},
        {"Pulay", 1},
        {"Broyden", 2}};

static std::unordered_map<std::string, int> charge_analysis = {
        {"None", 0},
        {"Voronoi", 1}};

static std::unordered_map<std::string, int> vdwdf_grid_type = {
        {"Coarse", 0},
        {"Fine", 1}};

static std::unordered_map<std::string, int> relax_mass = {
        {"Atomic", 0},
        {"Equal", 1}};

static std::unordered_map<std::string, int> md_temperature_control = {
        {"Nose Hoover Chains", 0},
        {"Anderson Rescaling", 1}};

static std::unordered_map<std::string, int> atomic_coordinate_type = {
        {"Cell Relative", 0},
        {"Absolute", 1}};

static std::unordered_map<std::string, int> calculation_mode = {
        {"Quench Electrons", 0},
        {"Relax Structure", 1},
        {"Constant Volume And Energy", 2},
        {"Constant Temperature And Energy", 3},
        {"Constant Pressure And Energy", 4},
        {"Plot", 5},
        {"Psi Plot", 6},
        {"Band Structure Only", 7},
        {"NEB Relax", 8},
        {"Dimer Relax", 9},
        {"TDDFT", 10},
        {"Exx Only", 11},
        {"STM", 13},
        {"NSCF", 14}
        };

static std::unordered_map<std::string, int> dos_method = {
        {"tetrahedra", 0},
        {"Gaussian", 1},
        };

static std::unordered_map<std::string, int> tetra_method = {
        {"Bloechl", 0},
        {"Linear", 1},
        {"Optimized", 2}
        };

static std::unordered_map<std::string, int> occupations_type = {
        {"Fixed", 0},
        {"Fermi Dirac", OCC_FD},
        {"Gaussian", OCC_GS},
        {"Error Function", OCC_EF},
        {"Cold Smearing", OCC_MV},
        {"MethfesselPaxton", OCC_MP},
        {"Tetrahedron", OCC_TETRA}
        };

static std::unordered_map<std::string, int> subdiag_driver = {
        {"lapack", SUBDIAG_LAPACK},
        {"scalapack", SUBDIAG_SCALAPACK},
        {"magma", SUBDIAG_MAGMA},
        {"cusolver", SUBDIAG_CUSOLVER},
        {"elpa", SUBDIAG_ELPA},
        {"rocsolver", SUBDIAG_ROCSOLVER},
        {"auto", SUBDIAG_AUTO}};

static std::unordered_map<std::string, int> kohn_sham_solver = {
        {"multigrid", MULTIGRID_SOLVER},
        {"davidson", DAVIDSON_SOLVER}};

static std::unordered_map<std::string, int> force_derivate_type = {
        {"wavefunction", WAVEFUNCTION_DERIVATIVE},
        {"projector", PROJECTOR_DERIVATIVE}};

static std::unordered_map<std::string, int> poisson_solver = {
        {"multigrid", MULTIGRID_SOLVER},
        {"pfft", POISSON_PFFT_SOLVER}};

static std::unordered_map<std::string, int> kpoint_units = {
        {"Reciprocal lattice", KPOINT_LATT_UNIT},
        {"2pi/alat", KPOINT_2pi_alat}};
 
static std::unordered_map<std::string, int> md_integration_order = {
        {"2nd Velocity Verlet", 0},
        {"3rd Beeman-Velocity Verlet", 1},
        {"5th Beeman-Velocity Verlet", 2}};

static std::unordered_map<std::string, int> interpolation_type = {
        {"Cubic Polynomial", CUBIC_POLYNOMIAL_INTERPOLATION},
        {"prolong", PROLONG_INTERPOLATION},
        {"FFT", FFT_INTERPOLATION}};

static std::unordered_map<std::string, int> start_mode = {
        {"Random Start", RANDOM_START},
        {"Restart From File", RESTART},
        {"LCAO Start", LCAO_START},
        {"FIREBALL Start", INIT_FIREBALL},
        {"Gaussian Start", INIT_GAUSSIAN},
        {"Start TDDFT", Start_TDDFT},
        {"Restart TDDFT", Restart_TDDFT},
        {"Modified LCAO Start", MODIFIED_LCAO_START}
        };

static std::unordered_map<std::string, int> ldaU_mode = {
        {"None", LDA_PLUS_U_NONE},
        {"Simple", LDA_PLUS_U_SIMPLE}};

static std::unordered_map<std::string, int> boundary_condition_type = {
        {"Periodic", 0}};
//        {"Cluster", 1},
//        {"Surface", 2}};

static std::unordered_map<std::string, int> z_average_output_mode = {
        {"None", 0},
        {"potential and charge density", 1},
        {"wave functions", 2}};

// The integer value maps to the position in reordered_xc_type
// immediately following
static std::unordered_map<std::string, int> exchange_correlation_type = {
        {"pz", 0},
        {"PZ", 0},
        {"LDA", 0},
        {"bp", 1},
        {"BP", 1},
        {"GGA XB CP", 1},
        {"pw91", 2},
        {"PW91", 2},
        {"GGA XP CP", 2},
        {"blyp", 3},
        {"BLYP", 3},
        {"GGA BLYP", 3},
        {"pbe", 4},
        {"PBE", 4},
        {"GGA PBE", 4},
        {"revpbe", 5},
        {"REVPBE", 5},
        {"pw86pbe", 6},
        {"PW86PBE", 6},
        {"b86bpbe", 7},
        {"pbesol", 8},
        {"PBESOL", 8},
        {"q2d", 9},
        {"Q2D", 9},
        {"hcth", 10},
        {"HCTH", 10},
        {"olyp", 11},
        {"wc", 12},
        {"sogga", 13},
        {"optbk88", 14},
        {"optb86b", 15},
        {"ev93", 16},
        {"tpss", 17},
        {"m06l", 18},
        {"tb09", 19},
        {"TB09", 19},
        {"MGGA TB09", 19},
        {"mgga tb09", 19},
        {"pbe0", 20},
        {"PBE0", 20},
        {"hse", 21},
        {"HSE", 21},
        {"b3lyp", 22},
        {"B3LYP", 22},
        {"gaupbe", 23},
        {"vdw-df", 24},
        {"VDW-DF", 24},
        {"vdw-df-cx", 25},
        {"VDW-DF-CX", 25},
        {"sla+pw+pbe+vdw1", 26},
        {"vdw-df-c09", 27},
        {"hartree-fock", 28},
        {"AUTO_XC", 99}};

// Internal types accepted by QE routines. Must map to preceding array
static std::string reordered_xc_type[] = {
        {"pz"},
        {"bp"},
        {"pw91"},
        {"blyp"},
        {"pbe"},
        {"revpbe"},
        {"pw86pbe"},
        {"b86bpbe"},
        {"pbesol"},
        {"q2d"},
        {"hcth"},
        {"olyp"},
        {"wc"},
        {"sogga"},
        {"optbk88"},
        {"optb86b"},
        {"ev93"},
        {"tpss"},
        {"m06l"},
        {"tb09"},
        {"pbe0"},
        {"hse"},
        {"b3lyp"},
        {"gaupbe"},
        {"vdw-df"},
        {"vdw-df-cx"},
        {"sla+pw+pbe+vdw1"},
        {"vdw-df-c09"},
        {"hf"}};

static std::unordered_map<std::string, int> relax_method = {
        {"Fast Relax", 0},
        {"FIRE", 1},
        {"Quick Min", 2},
        {"MD Min", 3},
        {"LBFGS", 4}};
static std::unordered_map<std::string, int> mg_method = {
        {"Steepest Descent", 0},
        {"Pulay", 1},
        {"Kain", 2}
        };
static std::unordered_map<std::string, int> orbital_layout = {
        {"Distribute", 0},
        {"Projection", 1}
        };

static std::unordered_map<std::string, int> energy_point_insert_mode = {
        {"None", 0},
        {"Simpson", 1},
        {"Sharp Peaks", 2}
        };

static std::unordered_map<std::string, int> exx_mode = {
        {"Distributed fft", EXX_DIST_FFT},
        {"Local fft", EXX_LOCAL_FFT},
        };

static std::unordered_map<std::string, int> exxdiv_treatment = {
        {"gygi-baldereschi", EXX_DIV_GYGI_BALDERESCHI},
        {"none", EXX_DIV_NONE},
        };
#endif
