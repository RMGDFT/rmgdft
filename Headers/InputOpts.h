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
        {"Hexagonal Rhombohedral (Trigonal)", 5},
        {"Tetragonal Primitive", 6},
        {"Tetragonal Body Centered", 7},
        {"Orthorhombic Primitive", 8},
        {"Orthorhombic Base Centered", 9},
        {"Orthorhombic Body Centered", 10},
        {"Orthorhombic Face Centered", 11},
        {"Monoclinic Primitive", 12},
        {"Monoclinic Base Centered", 13},
        {"Triclinic Primitive", 14}};

static std::unordered_map<std::string, int> energy_output_units = {
        {"Hartrees", 0},
        {"Rydbergs", 1}};

static std::unordered_map<std::string, int> discretization_type = {
        {"Mehrstellen", 0},
        {"Central", 1}};

static std::unordered_map<std::string, int> crds_units = {
        {"Bohr", 0},
        {"Angstrom", 1}};

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
        {"TDDFT", 10}
        };

static std::unordered_map<std::string, int> occupations_type = {
        {"Fixed", 0},
        {"Fermi Dirac", 1},
        {"Gaussian", 2},
        {"Error Function", 3},
        {"Cold Smearing", 4},
        {"MethfesselPaxton", 5}
        };

static std::unordered_map<std::string, int> subdiag_driver = {
        {"lapack", 0},
        {"scalapack", 1},
        {"magma", 2},
        {"auto", 3},
        {"elpa", 6}};

static std::unordered_map<std::string, int> kohn_sham_solver = {
        {"multigrid", MULTIGRID_SOLVER},
        {"davidson", DAVIDSON_SOLVER}};

static std::unordered_map<std::string, int> force_derivate_type = {
        {"wavefunction", WAVEFUNCTION_DERIVATIVE},
        {"projector", PROJECTOR_DERIVATIVE}};

static std::unordered_map<std::string, int> poisson_solver = {
        {"multigrid", MULTIGRID_SOLVER},
        {"pfft", POISSON_PFFT_SOLVER}};

static std::unordered_map<std::string, int> md_integration_order = {
        {"2nd Velocity Verlet", 0},
        {"3rd Beeman-Velocity Verlet", 1},
        {"5th Beeman-Velocity Verlet", 2}};

static std::unordered_map<std::string, int> interpolation_type = {
        {"Cubic Polynomial", CUBIC_POLYNOMIAL_INTERPOLATION},
        {"B-spline", BSPLINE_INTERPOLATION},
        {"prolong", PROLONG_INTERPOLATION},
        {"FFT", FFT_INTERPOLATION}};

static std::unordered_map<std::string, int> start_mode = {
        {"Random Start", 0},
        {"Restart From File", 1},
        {"LCAO Start", 2},
        {"FIREBALL Start", 3},
        {"Gaussian Start", 4},
        {"Start TDDFT", 5},
        {"Restart TDDFT", 6},
        {"Modified LCAO Start", 7}
        
        };

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
        {"vdw-df-c09"}};

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

static std::unordered_map<std::string, int> energy_point_insert_mode = {
        {"None", 0},
        {"Simpson", 1},
        {"Sharp Peaks", 2}
        };

#endif
