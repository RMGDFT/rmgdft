#ifndef RMG_InputOpts_H
#define RMG_InputOpts_H 1

#include <unordered_map>

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

static std::unordered_map<std::string, int> discretization_type = {
        {"Mehrstellen", 0},
        {"Central", 1}};

static std::unordered_map<std::string, int> crds_units = {
        {"Bohr", 0},
        {"Angstrom", 1}};

static std::unordered_map<std::string, int> charge_mixing_type = {
        {"Linear", 0},
        {"Pulay", 1}};

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
        {"Dimer Relax", 9}};

static std::unordered_map<std::string, int> occupations_type = {
        {"Fixed", 0},
        {"Fermi Dirac", 1},
        {"Gaussian", 2},
        {"Error Function", 3}};

static std::unordered_map<std::string, int> subdiag_driver = {
        {"lapack", 0},
        {"scalapack", 1},
        {"magma", 2}};

static std::unordered_map<std::string, int> md_integration_order = {
        {"2nd Velocity Verlet", 0},
        {"3rd Beeman-Velocity Verlet", 1},
        {"5th Beeman-Velocity Verlet", 2}};

static std::unordered_map<std::string, int> interpolation_type = {
        {"Cubic Polynomial", 0},
        {"B-spline", 1},
        {"prolong", 2}};

static std::unordered_map<std::string, int> start_mode = {
        {"Random Start", 0},
        {"Restart From File", 1},
        {"LCAO Start", 2},
        {"FIREBALL Start", 2},
        {"Gaussian Start", 3},
        {"Restart TDDFT", 4}
        
        };

static std::unordered_map<std::string, int> boundary_condition_type = {
        {"Periodic", 0}};
//        {"Cluster", 1},
//        {"Surface", 2}};

static std::unordered_map<std::string, int> z_average_output_mode = {
        {"None", 0},
        {"potential and charge density", 1},
        {"wave functions", 2}};

static std::unordered_map<std::string, int> exchange_correlation_type = {
        {"LDA", 0},
        {"GGA BLYP", 1},
        {"GGA XB CP", 2},
        {"GGA XP CP", 3},
        {"GGA PBE", 4},
        {"MGGA TB09", 5}};

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

#endif
