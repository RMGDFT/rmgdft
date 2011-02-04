/************************** SVN Revision Information **************************
 **    $Id: inputs.h 723 2007-02-02 09:12:25Z miro $    **
******************************************************************************/
 
/****f* QMD-MGDFT/options.h *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *                       Frisco Rose
 *                       Jerzy Bernholc
 * FUNCTION
 *   
 * INPUTS
 *
 * OUTPUT
 *  
 * PARENTS
 *
 * CHILDREN
 * 
 * SEE ALSO
 *  
 * SOURCE
 */
#include <stdio.h>
#include "md.h"

#define TRUE 1
#define FALSE 0

/*  generic name to status mapping for input options */
struct b_list bool_list[] =
{
    {
    "do_write_waves_to_file", FALSE},
    {
    "do_print_energy_and_eigenvalues", FALSE},
    {
    "do_initial_diagnolization", FALSE},
    {
    "do_md_randomize_velocity", FALSE},
    {
    "do_wait_flag", FALSE},
    {
    "do_sort_wavefunctions", FALSE},
    {
    "do_dx_movie", FALSE},
    {
    "do_rm_movie", FALSE},
    {
    "do_xbs_movie", FALSE},
    {
    "do_milliken_population", FALSE},
    {
    "do_restart_overwrite_occupation_numbers", FALSE},
    {
    "do_restart_coordinates_with_current_as_initial", FALSE},
    {
    "write_memory_report", FALSE},
    {
        "do_movable", FALSE},
    {
        "do_movable_orbital_centers", FALSE},
    {
        "do_spin_polarized", FALSE},
    {
        "metalic", TRUE},
    };

#define NBOOLS ( sizeof bool_list/ sizeof bool_list[0])
int nbools = NBOOLS;

/* List of input tags, all Flags initialized to FALSE except when 
   option is unique. */
struct t_list tag_list[] =
{
    {"start_mode", NULL, 3,  /* tagname, active option, option count */
        {
	    "Random Start", 
	    "Restart From File", 
	    "LCAO Start"}},
    {"boundary_condition_type", NULL, 3,
        {
	    "Periodic", 
	    "Cluster", 
	    "Surface"}},
    {"z_average_output_mode", NULL, 3,
	{
	    "None", 
	    "potential and charge density", 
	    "wave functions"}},
    {"exchange_correlation_type", NULL, 5,
        {
	    "LDA", 
	    "GGA BLYP", 
	    "GGA XB CP", 
	    "GGA XP CP", 
	    "GGA PBE"}},
    {"occupations_type", NULL, 4,
        {
	    "Fixed", 
	    "Fermi Dirac", 
	    "Gaussian", 
	    "Error Function"}},
    {"calculation_mode", NULL, 9,
        {
	    "Quench Electrons",
	    "Fast Relax",
	    "Constant Volume And Energy",
	    "Constant Temperature And Energy",
	    "Constant Pressure And Energy", 
	    "Constrained Fast Relax", 
        "Plot",
        "Psi Plot",
	    "Band Structure Only"}},
    {"md_temperature_control", NULL, 2,
        {
	    "Nose Hoover Chains", 
	    "Anderson Rescaling"}},
    {"md_integration_order", NULL, 3,
        {
	    "2nd Velocity Verlet", 
	    "3rd Beeman-Velocity Verlet", 
	    "5th Beeman-Velocity Verlet"}},
    {"bravais_lattice_type", NULL, 15,       /* Note, this is not cms compliant */
        {
	    "None",
	    "Cubic Primitive",
	    "Cubic Face Centered",
	    "Cubic Body Centered",
	    "Hexagonal Primitive",
	    "Hexagonal Rhombohedral (Trigonal)",
	    "Tetragonal Primitive",
	    "Tetragonal Body Centered",
	    "Orthorhombic Primitive",
	    "Orthorhombic Base Centered",
	    "Orthorhombic Body Centered",
	    "Orthorhombic Face Centered",
	    "Monoclinic Primitive", "Monoclinic Base Centered", "Triclinic Primitive"}},
    {"length_units", NULL, 2,      /* Note, this is not cms compliant */
	{
	    "Bohr", 
	    "Angstrom"}},
    {"atomic_coordinate_type", NULL, 2,      /* Note, this is not cms compliant */
	{
	    "Cell Relative", 
	    "Absolute"}},
    {"restart_coordinates_from", NULL, 2,    /* Note, this is not cms compliant */
	{
	    "Control File", 
	    "Wave File"}},
    {"interpolation_type", NULL, 2,  /* Note, this is not cms compliant */
	{
	    "Cubic Polynomial", 
	    "B-spline"}},

        {"mg_method", NULL, 3,
            {
                "Steepest Descent",
                "Pulay",
                "KAIN"}},

};

#define NTAGS ( sizeof tag_list/ sizeof tag_list[0])
int ntags = NTAGS;
