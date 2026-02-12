#pragma once
#ifndef RMG_InputOpts_H
#define RMG_InputOpts_H 1

#include <unordered_map>
#include <const.h>

extern std::unordered_map<std::string, int> bravais_lattice_type;
extern std::unordered_map<std::string, int> tddft_mode;
extern std::unordered_map<std::string, int> atomic_orbital_type;
extern std::unordered_map<std::string, int> internal_pseudo_type;
extern std::unordered_map<std::string, int> energy_output_units;
extern std::unordered_map<std::string, int> drho_precond_type;
extern std::unordered_map<std::string, int> crds_units;
extern std::unordered_map<std::string, int> vdw_corr;
extern std::unordered_map<std::string, int> lattice_units;
extern std::unordered_map<std::string, int> charge_mixing_type;
extern std::unordered_map<std::string, int> charge_analysis;
extern std::unordered_map<std::string, int> vdwdf_grid_type;
extern std::unordered_map<std::string, int> relax_mass;
extern std::unordered_map<std::string, int> md_temperature_control;
extern std::unordered_map<std::string, int> atomic_coordinate_type;
extern std::unordered_map<std::string, int> calculation_mode;
extern std::unordered_map<std::string, int> dos_method;
extern std::unordered_map<std::string, int> tetra_method;
extern std::unordered_map<std::string, int> occupations_type;
extern std::unordered_map<std::string, int> subdiag_driver;
extern std::unordered_map<std::string, int> kohn_sham_solver;
extern std::unordered_map<std::string, int> force_derivate_type;
extern std::unordered_map<std::string, int> poisson_solver;
extern std::unordered_map<std::string, int> kpoint_units;
extern std::unordered_map<std::string, int> md_integration_order;
extern std::unordered_map<std::string, int> interpolation_type;
extern std::unordered_map<std::string, int> start_mode;
extern std::unordered_map<std::string, int> ldaU_mode;
extern std::unordered_map<std::string, int> boundary_condition_type;
extern std::unordered_map<std::string, int> z_average_output_mode;
extern std::unordered_map<std::string, int> exchange_correlation_type;
extern std::string reordered_xc_type[];
extern std::unordered_map<std::string, int> relax_method;
extern std::unordered_map<std::string, int> mg_method;
extern std::unordered_map<std::string, int> orbital_layout;
extern std::unordered_map<std::string, int> energy_point_insert_mode;
extern std::unordered_map<std::string, int> exx_mode;
extern std::unordered_map<std::string, int> exxdiv_treatment;

#endif
