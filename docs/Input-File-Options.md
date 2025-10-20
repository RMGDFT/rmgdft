 -- A Real Space Multigrid Electronic structure code --
 --      More information at www.rmgdft.org          --


[\[Introduction \] ](#introduction)
[\[Control options\] ](#control-options)
[\[Cell parameter options\] ](#cell-parameter-options)
[\[Pseudopotential related options\] ](#pseudopotential-related-options)
[\[Kohn Sham solver options\] ](#kohn-sham-solver-options)
[\[Exchange correlation options\] ](#exchange-correlation-options)
[\[Orbital occupation options\] ](#orbital-occupation-options)
[\[Charge density mixing options\] ](#charge-density-mixing-options)
[\[Relaxation and Molecular dynamics options\] ](#relaxation-and-molecular-dynamics-options)
[\[Diagonalization options\] ](#diagonalization-options)
[\[Performance related options\] ](#performance-related-options)
[\[LDAU options\] ](#ldau-options)
[\[TDDFT related options\] ](#tddft-related-options)
[\[Poisson solver options\] ](#poisson-solver-options)
[\[Testing options\] ](#testing-options)
[\[Miscellaneous options\] ](#miscellaneous-options)


## Introduction

    The RMG input file consists of a set of key-value pairs of the form.
    
        name = "scalar"

    
    where scalar can be an integer, double or boolean value.

        period_of_diagonalization = "1"
        charge_density_mixing = "0.5"
        initial_diagonalization = "true"
    
    There are also strings and arrays which are delineated by double quotes so an
    integer array with three elements would be.
    
        processor_grid = "2 2 2"
    
    while a string example would be
    
        description = "64 atom diamond cell test run at the gamma point"
    
    strings can span multiple lines so the following would be valid as well.
    
        description = "64 atom diamond cell test run at gamma point
        using a Vanderbilt ultrasoft pseudopotential"
    
    string vectors span multiple lines and are used to enter items like a kpoint list one item per line.
    
        kpoints = "0.00  0.00  0.00   0.50
                   0.25  0.25  0.50   0.25
                   0.50  0.50  0.50   0.25"


Control options

    Key name:     AFM
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If true, anti-feromagnetic will be forced by symmetry operation if 
                  possible. 

    Key name:     BerryPhase
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  turn on/off Berry Phase calcualtion 

    Key name:     BerryPhaseCycle
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10
    Default:      1
    Description:  Berry Phase loop without updating rho and potentials 

    Key name:     BerryPhaseDirection
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2
    Default:      2
    Description:  Berry Phase direction: it will be efield direction when efield is 
                  non zero 

    Key name:     STM_bias
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "-1.0 1.0"
    Allowed:      
    Description:  Bias (in unit of Volt) for STM calculation 

    Key name:     STM_height
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "2.0 4.0"
    Allowed:      
    Description:  Height range for STM calculation 

    Key name:     a_length
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  First lattice constant. 

    Key name:     adaptive_cmix
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    10.000000
    Default:      0.000000e+00
    Description:  Manual setting for the adaptive interpolation parameter. 

    Key name:     afd_cfac
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    10.000000
    Default:      0.000000e+00
    Description:  Manual setting for the adaptive finite differencing parameter. 

    Key name:     b_length
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  Second lattice constant. 

    Key name:     c_length
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  Third lattice constant. 

    Key name:     calculation_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Quench Electrons"
    Allowed:      "NSCF" "Constant Temperature And Energy" "Constant Volume And 
                  Energy" "Dimer Relax" "TDDFT" "Relax Structure" "Constant Pressure 
                  And Energy" "Quench Electrons" "Plot" "Psi Plot" "Band Structure 
                  Only" "STM" "NEB Relax" "Exx Only" 
    Description:  Type of calculation to perform. 

    Key name:     cell_relax
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  flag to control unit cell relaxation 

    Key name:     coalesce_factor
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    16
    Default:      4
    Description:  Grid coalescing factor. 

    Key name:     coalesce_states
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to coalesce states. 

    Key name:     compressed_infile
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not parallel restart wavefunction file 
                  uses compressed format. 

    Key name:     compressed_outfile
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not parallel output wavefunction file 
                  uses compressed format. 

    Key name:     davidson_1stage_ortho
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag that improves davidson convergence for hard cases. Defaults 
                  to true but can be disabled for well behaved systems enabling 
                  higher performance. 

    Key name:     davidson_2stage_ortho
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Flag to futher improve davidson convergence for hard cases. 
                  Experimental. 

    Key name:     description
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Description of the run. 

    Key name:     drho_precond_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Resta"
    Allowed:      "Kerker" "Resta" 
    Description:  Density mixing preconditioner method. Resta or Kerker are 
                  supported. 

    Key name:     energy_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000e-20
    Max value:    1.000000e-07
    Default:      1.000000e-10
    Description:  The RMS value of the estimated change in the total energy per step 
                  where we assume self consistency has been achieved. 

    Key name:     energy_output_units
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Hartrees"
    Allowed:      "Rydbergs" "Hartrees" 
    Description:  Units to be used when writing energy values to the output file. 
                  Hartrees or Rydbergs are available. 

    Key name:     epsg_guard
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000e-05
    Default:      1.000000e-07
    Description:  GGA guard value for low density regions. 

    Key name:     exx_integrals_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "afqmc_rmg"
    Allowed:      
    Description:  File/path for exact exchange integrals. 

    Key name:     exx_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Local fft"
    Allowed:      "Local fft" "Distributed fft" 
    Description:  FFT mode for exact exchange computations. 

    Key name:     exxdiv_treatment
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "gygi-baldereschi"
    Allowed:      "none" "gygi-baldereschi" 
    Description:  Exact exchange method for handling exx divergence at G=0. 

    Key name:     freeze_ldaU_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      500
    Description:  freeze the ldaU occupations ns_occ after this step. 

    Key name:     gpu_managed_memory
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Some AMD and Nvidia GPUs support managed gou memory which is 
                  useful when GPU memory limits are exceeded. 

    Key name:     input_tddft_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave_tddft.out"
    Allowed:      
    Description:  Input file/path to read wavefunctions and other binary data from 
                  on a restart. 

    Key name:     input_wave_function_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave.out"
    Allowed:      
    Description:  Input file/path to read wavefunctions and other binary data from 
                  on a restart. 

    Key name:     interpolation_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "FFT"
    Allowed:      "FFT" "prolong" "Cubic Polynomial" 
    Description:  Interpolation method for transferring data between the potential 
                  grid and the wavefunction grid. Mostly for diagnostic purposes. 

    Key name:     lambda_max
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000
    Max value:    100.000000
    Default:      3.000000
    Description:  Chebyshev smoothing parameter. Don't change unless you know what 
                  you're doing. 

    Key name:     lambda_min
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.000000e+00
    Max value:    2.000000
    Default:      0.300000
    Description:  Chebyshev smoothing parameter. Don't change unless you know what 
                  you're doing. 

    Key name:     max_exx_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      100
    Description:  Maximum number of self consistent steps to perform with hybrid 
                  functionals. 

    Key name:     max_scf_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      100
    Description:  Maximum number of self consistent steps to perform. Inner loop for 
                  hybrid functionals. 

    Key name:     noncollinear
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If set true, noncollinear calculation. 

    Key name:     nvme_orbitals
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not orbitals should be mapped to disk. 

    Key name:     nvme_orbitals_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Orbitals/"
    Allowed:      
    Description:  File/path for runtime disk storage of orbitals. 

    Key name:     nvme_weights
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not projector weights should be mapped 
                  to disk. 

    Key name:     nvme_weights_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Weights/"
    Allowed:      
    Description:  File/path for disk storage of projector weights. 

    Key name:     nvme_work
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not work arrays should be mapped to 
                  disk. 

    Key name:     nvme_work_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Work/"
    Allowed:      
    Description:  File/path for disk storage of workspace. 

    Key name:     omp_threads_per_node
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    64
    Default:      0
    Description:  Number of Open MP threads each MPI process will use. A value of 0 
                  selects automatic setting. 

    Key name:     output_tddft_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave_tddft.out"
    Allowed:      
    Description:  Output file/path to store wavefunctions and other binary data. 

    Key name:     output_wave_function_file
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Waves/wave.out"
    Allowed:      
    Description:  Output file/path to store wavefunctions and other binary data. 

    Key name:     pseudo_dir
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "."
    Allowed:      
    Description:  Directory where pseudopotentials are stored. 

    Key name:     qfunction_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Qfunctions/"
    Allowed:      
    Description:  File/path for runtime disk storage of qfunctions. 

    Key name:     qmc_nband
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      0
    Description:  The number of band used in rmg-qmcpack interface. 

    Key name:     read_serial_restart
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Directs RMG to read from serial restart files. Normally used when 
                  changing the sprocessor topology used during a restart run 

    Key name:     rms_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000e-03
    Default:      1.000000e-07
    Description:  The RMS value of the change in the total potential from step to 
                  step where we assume self consistency has been achieved. 

    Key name:     semilocal_projectors
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    6
    Max value:    40
    Default:      10
    Description:  Controls the number of semilocal projectors. 

    Key name:     spinorbit
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If set true, spinorbit coupling calculation. 

    Key name:     start_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "LCAO Start"
    Allowed:      "Modified LCAO Start" "Restart TDDFT" "Start TDDFT" "FIREBALL 
                  Start" "Gaussian Start" "LCAO Start" "Restart From File" "Random 
                  Start" 
    Description:  Type of run. 

    Key name:     stress
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  flag to control stress cacluation 

    Key name:     stress_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    50.000000
    Default:      0.500000
    Description:  The stress criteria 

    Key name:     system_charge
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -unlimited
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  Number of excess holes in the system (useful for doped systems). 
                  Example: 2 means system is missing two electrons 

    Key name:     time_reversal
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if false, no k -> -k symmetry 

    Key name:     use_energy_correction
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Experimental energy correction term 

    Key name:     vdw_corr
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "None"
    Allowed:      "DFT-D3" "Grimme-D2" "DFT-D2" "None" 
    Description:  Type of vdw correction 

    Key name:     vdwdf_grid_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Fine"
    Allowed:      "Fine" "Coarse" 
    Description:  Type of grid to use when computing vdw-df correlation. 

    Key name:     vdwdf_kernel_filepath
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "vdW_kernel_table"
    Allowed:      
    Description:  File/path for vdW_kernel_table data. 

    Key name:     verbose
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag for writing out extra information 

    Key name:     wannier90
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  set up informations for wannier90 interface 

    Key name:     wannier90_scdm_mu
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -unlimited
    Max value:    unlimited
    Default:      0.000000e+00
    Description:  when wannier90 is used to build wannier functions, the energy 
                  window parameter 

    Key name:     write_data_period
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -5
    Max value:    50000
    Default:      -1
    Description:  How often to write checkpoint files during the initial quench in 
                  units of SCF steps. During structural relaxations of molecular 
                  dynamics checkpoints are written each ionic step. 

    Key name:     write_qmcpack_restart
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If true then a QMCPACK restart file is written as well as a serial 
                  restart file. 

    Key name:     write_qmcpack_restart_localized
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If true then a QMCPACK restart file for localized orbitals 

    Key name:     write_serial_restart
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  RMG normally writes parallel restart files. These require that 
                  restarts have the same processor topology. If write_serial_restart 
                  = "true" then RMG will also write a serial restart file that can 
                  be used with a different processor topology 

Cell parameter options

    Key name:     atomic_coordinate_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Absolute"
    Allowed:      "Absolute" "Cell Relative" 
    Description:  Flag indicated whether or not atomic coordinates are absolute or 
                  cell relative. 

    Key name:     bravais_lattice_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "None"
    Allowed:      "Orthorhombic Primitive" "Monoclinic Primitive" "Tetragonal 
                  Primitive" "Hexagonal Primitive" "Cubic Body Centered" "Triclinic 
                  Primitive" "Cubic Face Centered" "Cubic Primitive" "None" 
    Description:  Bravais Lattice Type. 

    Key name:     cell_movable
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 0 0 0 0 0 0 "
    Description:  9 numbers to control cell relaxation 

    Key name:     crds_units
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Bohr"
    Allowed:      "Angstrom" "Bohr" 
    Description:  Units for the atomic coordinates. 

    Key name:     frac_symmetry
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  For supercell calculation, one can disable the fractional 
                  translation symmetry 

    Key name:     grid_spacing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.350000
    Description:  Approximate grid spacing (bohr). 

    Key name:     kpoint_distribution
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -2147483647
    Max value:    2147483647
    Default:      -1
    Description:  This option affects kpoint parallelization. If there are M MPI 
                  procs then N = M/kpoint_distribution procs are assigned to each 
                  kpoint. M must be evenly divisible by kpoint_distribution. 

    Key name:     kpoint_is_shift
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 "
    Description:  Three-D layout of the kpoint shift. 

    Key name:     kpoint_mesh
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "1 1 1 "
    Description:  Three-D layout of the kpoint mesh. 

    Key name:     kpoint_units
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Reciprocal lattice"
    Allowed:      "2pi/alat" "Reciprocal lattice" 
    Description:  kpoint units for reading kpoint 

    Key name:     lattice_units
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Bohr"
    Allowed:      "Alat" "Angstrom" "Bohr" 
    Description:  Units for the lattice vectors 

    Key name:     lattice_vector
    Required:     no
    Key type:     double array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 0 0 0 0 0 0 "
    Description:  The simulation cell may be specified using either lattice vectors, 
                  a0, a1, a2 or by lattice constants and a bravais lattice type. If 
                  lattice vectors are used they should be entered as a 3x3 matrix. 

    Key name:     ldos_end_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "-1 -1 -1 "
    Description:  a ending grid point for ldos caclualtion 

    Key name:     ldos_start_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "-1 -1 -1 "
    Description:  a grid point for starting ldos caclualtion 

    Key name:     potential_grid_refinement
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    4
    Default:      2
    Description:  Ratio of the potential grid density to the wavefunction grid 
                  density. For example if the wavefunction grid is (72,72,72) and 
                  potential_grid_refinement = "2" then the potential grid would be 
                  (144,144,144). The default value is 2 but it may sometimes be 
                  beneficial to adjust this. (For USPP the minimum value is also 2 
                  and it cannot be set lower. NCPP can be set to 1). 

    Key name:     processor_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "1 1 1 "
    Description:  Three-D (x,y,z) layout of the MPI processes. 

    Key name:     sts_end_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "-1 -1 -1 "
    Description:  a ending grid point for sts caclualtion 

    Key name:     sts_start_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "-1 -1 -1 "
    Description:  a grid point for starting sts caclualtion 

    Key name:     use_symmetry
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2
    Default:      2
    Description:  0: never use symmetry, 1: always use symmetry, 

    Key name:     wavefunction_grid
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "1 1 1 "
    Description:  Three-D (x,y,z) dimensions of the grid the wavefunctions are 
                  defined on. 

Pseudopotential related options

    Key name:     all_electron_parm
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: Yes
    Min value:    1
    Max value:    12
    Default:      4
    Description:  Gygi all electron parameter. 

    Key name:     atomic_orbital_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "delocalized"
    Allowed:      "delocalized" "localized" 
    Description:  Atomic Orbital Type. Choices are localized and delocalized. 

    Key name:     energy_cutoff_parameter
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.600000
    Max value:    1.000000
    Default:      0.800000
    Description:  

    Key name:     filter_factor
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.060000
    Max value:    1.000000
    Default:      1.000000
    Description:  Filtering factor. 

    Key name:     internal_pseudo_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "sg15"
    Allowed:      "all_electron" "nc_standard" "nc_accuracy" "sg15" "ultrasoft" 
    Description:  Internal pseudopotential type. Choices are sg15, ultrasoft, 
                  nc_accuracy or all_electron 

    Key name:     localize_localpp
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  The local potential associated with a particular ion also decays 
                  rapidly in real-space with increasing r. As with beta projectors 
                  truncating the real-space representation for large cells can lead 
                  to significant computational savings with a small loss of accuracy 
                  but it should be set to false for small cells. 

    Key name:     localize_projectors
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  The Beta function projectors for a particular ion decay rapidly in 
                  real-space with increasing r. For large cells truncating the 
                  real-space representation of the projector can lead to significant 
                  computational savings with a small loss of accuracy. For smaller 
                  cells the computational cost is the same for localized or 
                  delocalized projectors so it is better to set localize_projectors 
                  to false. 

    Key name:     max_nlradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    2.000000
    Max value:    10000.000000
    Default:      10000.000000
    Description:  maximum radius for non-local projectors 

    Key name:     max_qradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    2.000000
    Max value:    10000.000000
    Default:      10000.000000
    Description:  maximum radius for qfunc in ultra-pseudopotential 

    Key name:     min_nlradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000
    Max value:    10000.000000
    Default:      2.000000
    Description:  minimum radius for non-local projectors 

    Key name:     min_qradius
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000
    Max value:    10000.000000
    Default:      2.000000
    Description:  minimum radius for qfunc in ultra-pseudopotential 

    Key name:     projector_expansion_factor
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.500000
    Max value:    3.000000
    Default:      1.000000
    Description:  When using localized projectors the radius can be adjusted with 
                  this parameter. 

    Key name:     pseudopotential
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  External pseudopotentials may be specfied with this input key. The 
                  format uses the atomic symbol followed by the pseudopotential file 
                  name. pseudopotential = "Ni Ni.UPF O O.UPF" 

    Key name:     use_bessel_projectors
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  When a semi-local pseudopotential is being used projectors will be 
                  generated using Bloechl's procedure with Bessel functions as the 
                  basis set if this is true. 

Kohn Sham solver options

    Key name:     davidson_max_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    20
    Default:      8
    Description:  Maximum number of iterations for davidson diagonalization. 

    Key name:     davidson_multiplier
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    6
    Default:      0
    Description:  The davidson solver expands the eigenspace with the maximum 
                  expansion factor being set by the value of davidson_multiplier. 
                  Larger values often lead to faster convergence but because the 
                  computational cost of the davidson diagonalization step scales as 
                  the cube of the number of eigenvectors the optimal value based on 
                  the fastest time to solution depends on the number of orbitals. If 
                  not specified explicitly or set to 0 RMG uses the following 
                  algorithm to set the value. 
                  
                  Number of orbitals <= 600 davidson_multiplier= "4" 
                  600 < Number of orbitals <= 900 davidson_multiplier = "3" 
                  Number of orbitals > 900 davidson_multiplier = "2" 
                  
                  For very large problems the N^3 scaling makes even a factor of 2 
                  prohibitively costly and the multigrid solver is a better choice. 

    Key name:     davidson_premg
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    8
    Default:      0
    Description:  If the davidson solver is selected this parameter controls the 
                  number of multigrid steps to use before enabling davidson. 

    Key name:     kohn_sham_coarse_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.200000
    Default:      1.000000
    Description:  Time step to use in the kohn-sham multigrid solver on the coarse 
                  levels. 

    Key name:     kohn_sham_fd_order
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    6
    Max value:    12
    Default:      8
    Description:  RMG uses finite differencing to represent the kinetic energy 
                  operator and the accuracy of the representation is controllable by 
                  the kohn_sham_fd_order parameter. The default is 8 and is fine for 
                  most purposes but higher accuracy is obtainable with 10th or 12th 
                  order at the cost of some additional computational expense. 

    Key name:     kohn_sham_mg_levels
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -1
    Max value:    6
    Default:      -1
    Description:  Number of multigrid levels to use in the kohn-sham multigrid 
                  preconditioner. 

    Key name:     kohn_sham_mg_timestep
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.000000e+00
    Max value:    2.000000
    Default:      1.000000
    Description:  timestep for multigrid correction. 

    Key name:     kohn_sham_mucycles
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    6
    Default:      3
    Description:  Number of mu (also known as W) cycles to use in the kohn-sham 
                  multigrid preconditioner. 

    Key name:     kohn_sham_post_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    0
    Max value:    5
    Default:      1
    Description:  Number of global grid post-smoothing steps to perform after a 
                  multigrid preconditioner iteration. 

    Key name:     kohn_sham_pre_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    5
    Default:      2
    Description:  Number of global grid pre-smoothing steps to perform before a 
                  multigrid preconditioner iteration. 

    Key name:     kohn_sham_solver
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "davidson"
    Allowed:      "davidson" "multigrid" 
    Description:  RMG supports a pure multigrid Kohn-Sham solver as well as a 
                  multigrid preconditioned davidson solver. The davidson solver is 
                  usually better for smaller problems with the pure multigrid solver 
                  often being a better choice for very large problems. 

    Key name:     kohn_sham_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.000000e+00
    Max value:    2.000000
    Default:      0.660000
    Description:  Smoothing timestep to use on the fine grid in the the kohn-sham 
                  multigrid preconditioner. 

    Key name:     prolong_order
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    0
    Max value:    12
    Default:      10
    Description:  Debug option that controls interpolation order used to form the 
                  charge density and to compute the kinetic component of stress. If 
                  a value of 0 is selected then an FFT will be used. 

    Key name:     unoccupied_tol_factor
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000
    Max value:    100000.000000
    Default:      1000.000000
    Description:  When using the Davidson Kohn-Sham solver unoccupied states are 
                  converged to a less stringent tolerance than occupied orbitals 
                  with the ratio set by this parameter. 

    Key name:     use_block_diag
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Flag indicating whether or not to use block diagonalization. 

    Key name:     use_rmm_diis
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not to use the RMM-DIIS algorithm in 
                  the mulgrid solver. 

Exchange correlation options

    Key name:     exchange_correlation_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "AUTO_XC"
    Allowed:      "hartree-fock" "VDW-DF" "vdw-df" "gaupbe" "HSE" "PBE0" "mgga tb09" 
                  "TB09" "VDW-DF-CX" "optbk88" "AUTO_XC" "BLYP" "hse" "tb09" "GGA 
                  BLYP" "GGA XB CP" "ev93" "GGA XP CP" "Q2D" "pbe0" "PW91" "pbe" 
                  "bp" "vdw-df-c09" "optb86b" "pbesol" "pw91" "b86bpbe" "BP" "blyp" 
                  "PBE" "pz" "PZ" "LDA" "hcth" "revpbe" "m06l" "REVPBE" "vdw-df-cx" 
                  "pw86pbe" "B3LYP" "PW86PBE" "PBESOL" "sla+pw+pbe+vdw1" "sogga" 
                  "q2d" "MGGA TB09" "GGA PBE" "wc" "tpss" "HCTH" "b3lyp" "olyp" 
    Description:  Most pseudopotentials specify the exchange correlation type they 
                  were generated with and the default value of AUTO_XC means that 
                  the type specified in the pseudopotial is what RMG will use. That 
                  can be overridden by specifying a value here. 

    Key name:     exx_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000e-12
    Max value:    1.000000e-06
    Default:      1.000000e-09
    Description:  Convergence criterion for the EXX delta from step to step where we 
                  assume EXX consistency has been achieved. 

    Key name:     exx_fraction
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -1.000000
    Max value:    1.000000
    Default:      -1.000000
    Description:  when hybrid functional is used, the fraction of Exx 

    Key name:     vexx_fft_threshold
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000e-14
    Max value:    0.100000
    Default:      1.000000e-14
    Description:  The value for the EXX delta where we switch from single to double 
                  precision ffts. Single precision ffts are generally accurate 
                  enough. 

    Key name:     x_gamma_extrapolation
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set true, use exx extrapolation to gamma 

Orbital occupation options

    Key name:     MP_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    5
    Default:      2
    Description:  Order of Methefessel Paxton occupation. 

    Key name:     dos_broading
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.100000
    Description:  For DOS with Gaussian broading method 

    Key name:     dos_method
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "tetrahedra"
    Allowed:      "Gaussian" "tetrahedra" 
    Description:  tetrahedra or gauss smearing method for DOS calculation 

    Key name:     occupation_electron_temperature_eV
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    2.000000
    Default:      0.040000
    Description:  Target electron temperature when not using fixed occupations. 

    Key name:     occupation_number_mixing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      1.000000
    Description:  Mixing parameter for orbital occupations when not using fixed 
                  occupations. 

    Key name:     occupations_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Fermi Dirac"
    Allowed:      "Tetrahedron" "MethfesselPaxton" "Cold Smearing" "Error Function" 
                  "Gaussian" "Fermi Dirac" "Fixed" 
    Description:  RMG supports several different ways of specifying orbital 
                  occupations. For a spin polarized system one may specify the 
                  occupations for up and down separately. In the case of a non-zero 
                  electronic temperature these will be adjusted as the calculation 
                  proceeds based on this setting. 

    Key name:     states_count_and_occupation
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Occupation string for states. Format for a system with 240 
                  electrons and 20 unoccupied states would be. "120 2.0 20 0.0" 

    Key name:     states_count_and_occupation_spin_down
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Occupation string for spin down states. Format is the same as for 
                  states_count_and_occupation. Total number of states must match 
                  spin up occupation string. 

    Key name:     states_count_and_occupation_spin_up
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Occupation string for spin up states. Format is the same as for 
                  states_count_and_occupation. Total number of states must match 
                  spin down occupation string. 

    Key name:     tetra_method
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Bloechl"
    Allowed:      "Optimized" "Linear" "Bloechl" 
    Description:  tetrahedron method to use 

    Key name:     unoccupied_states_per_kpoint
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -2147483647
    Max value:    2147483647
    Default:      -1
    Description:  The number of unoccupied orbitals. A value that is 15-20% of the 
                  number of occupied orbitals generally works well. 

Charge density mixing options

    Key name:     charge_broyden_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    12
    Default:      10
    Description:  Number of previous steps to use when Broyden mixing is used to 
                  update the charge density. 

    Key name:     charge_density_mixing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.500000
    Description:  Proportion of the current charge density to replace with the new 
                  density after each scf step when linear mixing is used. 

    Key name:     charge_mixing_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Auto"
    Allowed:      "Broyden" "Auto" "Pulay" "Linear" 
    Description:  RMG supports Broyden, Pulay and Linear mixing When the davidson 
                  Kohn-Sham solver is selected Broyden or Pulay are preferred. For 
                  the multigrid solver Broyden is usually the best choice. 

    Key name:     charge_pulay_Gspace
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set true, charge density mixing the residual in G space 

    Key name:     charge_pulay_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10
    Default:      5
    Description:  Number of previous steps to use when Pulay mixing is used to 
                  update the charge density. 

    Key name:     charge_pulay_refresh
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      100
    Description:  charge Pulay mixing reset steps. 

    Key name:     drho_precond
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set true, charge density residual is preconded with 
                  q^2/(q^2+q0^2) 

    Key name:     drho_precond_q0
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    10.000000
    Default:      0.250000
    Description:  Kerker type preconditioning the charge density residual by 
                  q^2/(q^2+q0^2) See Kresse and Furthmueller, Computational 
                  Materials Science 6 (1996) 15-50 

    Key name:     ldau_mixing
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      1.000000
    Description:  Proportion of the current ldau occupation to replace with the new 
                  ones after each scf step when linear mixing is used. 

    Key name:     ldau_mixing_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Linear"
    Allowed:      "Broyden" "Auto" "Pulay" "Linear" 
    Description:  RMG supports Pulay and Linear mixing for DFT+U occupation mixing 

    Key name:     ldau_pulay_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10
    Default:      5
    Description:  Number of previous steps to use when Pulay mixing is used to 
                  update the ldau occupation . 

    Key name:     ldau_pulay_refresh
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      100
    Description:  ldau pulay mixing reset steps 

    Key name:     ldau_pulay_scale
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      1.000000
    Description:  

    Key name:     potential_acceleration_constant_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    4.000000
    Default:      1.000000
    Description:  When set to a non-zero value this parameter causes RMG to perform 
                  a band by band update of the self-consistent potential during the 
                  course of an SCF step when the multigrid kohn_sham_solver is 
                  chosen. This means that updates to the lower energy orbitals are 
                  incorporated into the SCF potential seen by the higher energy 
                  orbitals as soon as they are computed. This can lead to faster 
                  convergence and better stability for many systems. The option 
                  should only be used with Linear mixing. Even when the davidson 
                  solver is chosen this parameter may be used since the first few 
                  steps with davidson usually uses the multigrid solver. 

    Key name:     resta_beta
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1000.000000
    Default:      0.000000e+00
    Description:  Beta parameter for resta charge density preconditioning. The 
                  default value of 0.0 means the value should be autmatically 
                  determined. 

Relaxation and Molecular dynamics options

    Key name:     dynamic_time_counter
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     dynamic_time_delay
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    5
    Default:      5
    Description:  

    Key name:     force_grad_order
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    12
    Default:      8
    Description:  Atomic forces may be computed to varying degrees of accuracy 
                  depending on the requirements of a specific problem. A value of 0 
                  implies highest accuracy which is obtained by using FFTs in place 
                  of finite differencing. 

    Key name:     ionic_time_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      50.000000
    Description:  Ionic time step for use in molecular dynamics and structure 
                  optimizations. 

    Key name:     ionic_time_step_decrease
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.500000
    Description:  Factor by which ionic timestep is decreased when dynamic timesteps 
                  are enabled. 

    Key name:     ionic_time_step_increase
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000
    Max value:    3.000000
    Default:      1.100000
    Description:  Factor by which ionic timestep is increased when dynamic timesteps 
                  are enabled. 

    Key name:     max_ionic_time_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    150.000000
    Default:      150.000000
    Description:  Maximum ionic time step to use for molecular dynamics or 
                  structural optimizations. 

    Key name:     max_md_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      100
    Description:  Maximum number of molecular dynamics steps to perform. 

    Key name:     md_integration_order
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "5th Beeman-Velocity Verlet"
    Allowed:      "5th Beeman-Velocity Verlet" "3rd Beeman-Velocity Verlet" "2nd 
                  Velocity Verlet" 
    Description:  Integration order for molecular dynamics. 

    Key name:     md_nose_oscillation_frequency_THz
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      15.590000
    Description:  

    Key name:     md_number_of_nose_thermostats
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    5
    Default:      5
    Description:  Number of Nose thermostats to use during Constant Volume and 
                  Temperature MD. 

    Key name:     md_randomize_velocity
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  The initial ionic velocities for a molecular dyanamics run are 
                  randomly initialized to the target temperature. 

    Key name:     md_temperature
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      300.000000
    Description:  Target MD Temperature. 

    Key name:     md_temperature_control
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Nose Hoover Chains"
    Allowed:      "Anderson Rescaling" "Nose Hoover Chains" 
    Description:  Type of temperature control method to use in molecular dynamics. 

    Key name:     relax_dynamic_timestep
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to use dynamic timesteps in 
                  relaxation mode. 

    Key name:     relax_mass
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Atomic"
    Allowed:      "Equal" "Atomic" 
    Description:  Mass to use for structural relaxation, either atomic masses, or 
                  the mass of carbon for all atoms. 

    Key name:     relax_max_force
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      2.500000e-03
    Description:  Force value at which an ionic relaxation is considered to be 
                  converged. 

    Key name:     relax_method
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "LBFGS"
    Allowed:      "LBFGS" "MD Min" "Quick Min" "FIRE" "Fast Relax" 
    Description:  Type of relaxation method to use for structural optimizations. 

    Key name:     renormalize_forces
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not to renormalize forces. 

    Key name:     tddft_frequency
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.200000
    Description:  TDDFT frequency for use in TDDFT vector potential mode 

    Key name:     tddft_time_step
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      0.200000
    Description:  TDDFT time step for use in TDDFT mode 

Diagonalization options

    Key name:     extra_random_lcao_states
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      0
    Description:  LCAO (Linear Combination of Atomic Orbitals) is the default 
                  startup method for RMG. The atomic orbitals are obtained from the 
                  pseudpotentials but in some cases better convergence may be 
                  obtained by adding extra random wavefunctions in addition to the 
                  atomic orbitals. 

    Key name:     folded_spectrum
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  When the number of eigenvectors is large using folded_spectrum is 
                  substantially faster than standard diagonalization. It also tends 
                  to converge better for metallic systems. It works with the 
                  multigrid kohn_sham_solver but not the davidson solver. 

    Key name:     folded_spectrum_iterations
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    0
    Max value:    20
    Default:      2
    Description:  Number of folded spectrum iterations to perform. 

    Key name:     folded_spectrum_width
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.100000
    Max value:    1.000000
    Default:      0.300000
    Description:  Submatrix width to use as a fraction of the full spectrum. The 
                  folded spectrum width ranges from 0.10 to 1.0. For insulators and 
                  semiconductors a value of 0.3 is appropriate. For metals values 
                  between 0.15 to 0.2 tend to be better. The default value is 0.3 

    Key name:     initial_diagonalization
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Perform initial subspace diagonalization. 

    Key name:     period_of_diagonalization
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      1
    Description:  Diagonalization period (per scf step). Mainly for debugging and 
                  should not be changed for production. 

    Key name:     scalapack_block_factor
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    4
    Max value:    2147483647
    Default:      64
    Description:  Block size to use with scalapack. Optimal value is dependent on 
                  matrix size and system hardware. 

    Key name:     subdiag_driver
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "auto"
    Allowed:      "auto" "rocsolver" "elpa" "cusolver" "magma" "scalapack" "lapack" 
    Description:  Driver type used for subspace diagonalization of the eigenvectors. 

    Key name:     subdiag_groups
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    16
    Default:      1
    Description:  Number of scalapack or elpa groups. 

Performance related options

    Key name:     fd_allocation_limit
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1024
    Max value:    262144
    Default:      65536
    Description:  Allocation sizes in finite difference routines less than this 
                  value are stack rather than heap based. 

    Key name:     mpi_queue_mode
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Use mpi queue mode. 

    Key name:     non_local_block_size
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    64
    Max value:    40000
    Default:      512
    Description:  Block size to use when applying the non-local and S operators. A 
                  value at least as large as the number of wavefunctions produces 
                  better performance but requires more memory. 

    Key name:     preconditioner_threshold
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    1.000000e-09
    Max value:    0.100000
    Default:      0.100000
    Description:  The RMS value of the change in the total potential where we switch 
                  the preconditioner from single to double precision. 

    Key name:     require_huge_pages
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  If set RMG assumes that sufficient huge pages are available. Bad 
                  things may happen if this is not true. 

    Key name:     rmg_threads_per_node
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    0
    Max value:    64
    Default:      0
    Description:  Number of Multigrid/Davidson threads each MPI process will use. A 
                  value of 0 means set automatically. 

    Key name:     spin_manager_thread
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "true"
    Description:  When mpi_queue_mode is enabled the manager thread spins instead of 
                  sleeping. 

    Key name:     spin_worker_threads
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "true"
    Description:  When mpi_queue_mode is enabled the worker threads spin instead of 
                  sleeping. 

    Key name:     state_block_size
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      64
    Description:  State block size used in nlforce. Larger values require more 
                  memory but can 

    Key name:     use_alt_zgemm
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to use alternate zgemm 
                  implementation. 

    Key name:     use_async_allreduce
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  RMG uses MPI_Allreduce function calls in several places and for 
                  large problems these can account for a significant fraction of the 
                  total run time. In most cases using the asynchronous MPI versions 
                  of the functions is faster but this is not true for all platforms 
                  and in that casesetting this flag to false can improve 
                  performance. 

    Key name:     use_hwloc
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Use internal hwloc setup if available. If both this and use_numa 
                  are true hwloc takes precedence. 

    Key name:     use_numa
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  Numa stands for Non Uniform Memory Access and means that the main 
                  memory of a computer is organized into seperate distinct banks. 
                  Each bank is then attached to a CPU core or group of cores and 
                  while all cores can normally access all banks the access speed is 
                  faster for directly attached banks. Ensuring that individual CPU 
                  cores mostly access memory in banks they are directly attached to 
                  can have a large impact on performance. Process mapping that does 
                  this can normally be done when jobs are submitted and run via 
                  arguments to mpirun/mpiexec but if this is not done RMG will 
                  attempt to provide an optimal mapping if use_numa is set to true. 

LDAU options

    Key name:     Hubbard_U
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Hubbard U parameter for each atomic species using the format 
                  Hubbard_U="Ni 6.5 3d 0.0 0.0 0.0" 

    Key name:     ldaU_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "None"
    Allowed:      "Simple" "None" 
    Description:  Type of lda+u implementation. 

    Key name:     ldaU_radius
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000
    Max value:    12.000000
    Default:      9.000000
    Description:  Max radius of atomic orbitals to be used in LDA+U projectors. 

TDDFT related options

    Key name:     restart_tddft
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  restart TDDFT 

    Key name:     tddft_gpu
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  use gpu for ELYDYN or not 

    Key name:     tddft_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "electric field"
    Allowed:      "vector potential" "point charge" "electric field" 
    Description:  TDDFT mode 

    Key name:     tddft_noscf
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  TDDFT run read data directly from the last scf job 

    Key name:     tddft_qgau
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      1.000000
    Description:  Gaussian parameter for point charge to Gaussian charge 

    Key name:     tddft_qpos
    Required:     no
    Key type:     double array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 "
    Description:  cartesian coordinate of the point charge for tddft 

    Key name:     tddft_start_state
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      0
    Description:  the starting state to use in tddft dynamics 

    Key name:     tddft_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      2000
    Description:  Maximum number of tddft steps to perform. 

Poisson solver options

    Key name:     hartree_max_sweeps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    5
    Max value:    100
    Default:      10
    Description:  Maximum number of hartree iterations to perform per scf step. 

    Key name:     hartree_min_sweeps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    5
    Default:      5
    Description:  Minimum number of hartree iterations to perform per scf step. 

    Key name:     hartree_rms_ratio
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1000.000000
    Max value:    unlimited
    Default:      100000.000000
    Description:  Ratio between target RMS for get_vh and RMS total potential. 

    Key name:     poisson_coarse_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.400000
    Max value:    1.000000
    Default:      0.800000
    Description:  Time step to use in the poisson multigrid solver on the coarse 
                  levels. 

    Key name:     poisson_coarsest_steps
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    10
    Max value:    100
    Default:      25
    Description:  Number of smoothing steps to use on the coarsest level in the 
                  hartree multigrid solver. 

    Key name:     poisson_finest_time_step
    Required:     no
    Key type:     double
    Expert:       Yes
    Experimental: No
    Min value:    0.400000
    Max value:    1.000000
    Default:      1.000000
    Description:  Time step to use in the poisson multigrid solver on the finest 
                  level. 

    Key name:     poisson_mg_levels
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -1
    Max value:    6
    Default:      -1
    Description:  Number of multigrid levels to use in the hartree multigrid solver. 

    Key name:     poisson_mucycles
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    4
    Default:      3
    Description:  Number of mu (also known as W) cycles to use in the hartree 
                  multigrid solver. 

    Key name:     poisson_post_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    6
    Default:      1
    Description:  Number of global hartree grid post-smoothing steps to perform 
                  after a multigrid iteration. 

    Key name:     poisson_pre_smoothing
    Required:     no
    Key type:     integer
    Expert:       Yes
    Experimental: No
    Min value:    1
    Max value:    6
    Default:      2
    Description:  Number of global hartree grid pre-smoothing steps to perform 
                  before a multigrid iteration. 

    Key name:     poisson_solver
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "pfft"
    Allowed:      "pfft" "multigrid" 
    Description:  poisson solver. 



    Key name:     charge_analysis_period
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    500
    Default:      10
    Description:  How often to perform and write out charge analysis. 

    Key name:     output_rho_xsf
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Generate xsf format for electronic density. 

    Key name:     write_eigvals_period
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    100
    Default:      5
    Description:  How often to output eigenvalues in units of scf steps. 

    Key name:     write_orbital_overlaps
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If true the orbital overlap matrix from successive MD steps is 
                  written. 

    Key name:     write_pdos
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag to write partial density of states. 

    Key name:     write_pseudopotential_plots
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Flag to indicate whether or not to write pseudopotential plots. 

Testing options

    Key name:     test_bond_length
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    20.000000
    Default:      nan
    Description:  Expected dimer bond length for testing. 

    Key name:     test_bond_length_tolerance
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000e-04
    Max value:    0.100000
    Default:      1.000000e-03
    Description:  Test bond length tolerance. 

    Key name:     test_energy
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -1000000000.000000
    Max value:    1000000000.000000
    Default:      nan
    Description:  Expected final energy for testing. 

    Key name:     test_energy_tolerance
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    1.000000e-08
    Max value:    1.000000e-04
    Default:      1.000000e-07
    Description:  Test final energy tolerance. 

    Key name:     test_steps
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    1000
    Default:      0
    Description:  Expected number of scf steps for testing. 

    Key name:     test_steps_tolerance
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    10
    Default:      1
    Description:  Test scf steps tolerance. 

Miscellaneous options

    Key name:     E_POINTS
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    2147483647
    Default:      201
    Description:  number of E points for ldos and sts calculation 

    Key name:     Emax
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -100.000000
    Max value:    100.000000
    Default:      0.000000e+00
    Description:  

    Key name:     Emin
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    -100.000000
    Max value:    100.000000
    Default:      -6.000000
    Description:  

    Key name:     ExxCholMax
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    64
    Default:      8
    Description:  maximum number of Exx integral cholesky vectors 

    Key name:     ExxIntCholosky
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set true, Exx integrals are Cholesky factorized to 3-index 

    Key name:     adaptive_convergence
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Parameters that control initial SCF convergence are adaptively 
                  modified. Affected parameters include density mixing and 
                  preconditioning and electron temperature. 

    Key name:     alt_laplacian
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "true"
    Description:  Flag indicating whether or not to use alternate laplacian weights 
                  for some operators. 

    Key name:     boundary_condition_type
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Periodic"
    Allowed:      "Periodic" 
    Description:  Boundary condition type Only periodic is currently implemented. 

    Key name:     charge_analysis
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "Voronoi"
    Allowed:      "Voronoi" "None" 
    Description:  Type of charge analysis to use. Only Voronoi deformation density 
                  is currently available. 

    Key name:     cube_pot
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, total potential is printed out in cube format 

    Key name:     cube_rho
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, charge density is printed out in cube format 

    Key name:     cube_states_list
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  plot the states listed here 

    Key name:     cube_vh
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  if set true, hatree potential is printed out in cube format 

    Key name:     dftd3_version
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    2
    Max value:    6
    Default:      3
    Description:  Grimme's DFT-D3 versions, 

    Key name:     dipole_correction
    Required:     no
    Key type:     integer array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 "
    Description:  (1,1,1) for molecule, dipole correction in all directions. 

    Key name:     dipole_moment
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Turns on calculation of dipole moment for the entire cell. 

    Key name:     ecutrho
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    10000.000000
    Default:      0.000000e+00
    Description:  ecut for rho in unit of Ry. 

    Key name:     ecutwfc
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    10000.000000
    Default:      0.000000e+00
    Description:  ecut for wavefunctions in unit of Ry. 

    Key name:     electric_field
    Required:     no
    Key type:     double array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 "
    Description:  Components of the electric field in reciprocal lattice direction 
                  and unit of Ha/bohr. 

    Key name:     electric_field_tddft
    Required:     no
    Key type:     double array
    Expert:       No
    Experimental: No
    Default:      "0 0 0 "
    Description:  the electric field for TDDFT in reciprocal lattice direction and 
                  unit of Ha/bohr 

    Key name:     equal_initial_density
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Specifies whether to set initial up and down density to be equal. 

    Key name:     exx_int_flag
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  If set true, calculate the exact exchange integrals. 

    Key name:     fast_density
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Use a faster but less accurate method to generate the charge 
                  density from the electronic wavefunctions. As the cutoff 
                  (grid-density) increases this method improves in accuracy. This 
                  option should be set to false if you receive warnings about 
                  negative charge densities after interpolation. 

    Key name:     freeze_occupied
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Flag indicating whether or not to freeze the density and occupied 
                  orbitals after a restart. 

    Key name:     gw_residual_convergence_criterion
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: Yes
    Min value:    1.000000e-14
    Max value:    4.000000e-04
    Default:      1.000000e-06
    Description:  The max value of the residual for unoccupied orbitals when 
                  performing a GW calculation. 

    Key name:     gw_residual_fraction
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: Yes
    Min value:    0.000000e+00
    Max value:    1.000000
    Default:      0.900000
    Description:  The residual value specified by gw_residual_convergence_criterion 
                  is applied to this fraction of the total spectrum. 

    Key name:     kohn_sham_ke_fft
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Special purpose flag which will force use of an FFT for the 
                  kinetic energy operator. 

    Key name:     kpoints
    Required:     yes
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  Normally kpoints are specified using the kpoint_mesh and 
                  kpoint_is_shift options but one can also enter a list of kpoints 
                  and their weights with this option. If kpoint_mesh is not 
                  specified or this is a bandstructure calculation this is required 
                  otherwise it is optional. 

    Key name:     kpoints_bandstructure
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      ""
    Allowed:      
    Description:  List of kpoints to use in a bandstructure calculation. For more 
                  detailed information look at the github wiki page on kpoint 
                  calculations. 

    Key name:     lcao_use_empty_orbitals
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "false"
    Description:  Some pseudopotentials contain unbound atomic orbitals and this 
                  flag indicates whether or not they should be used for LCAO starts. 

    Key name:     md_steps_offset
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     num_wanniers
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    2147483647
    Default:      0
    Description:  number of wannier functions to be used in wannier90 

    Key name:     rmg2bgw
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Write wavefunction in G-space to BerkeleyGW WFN file. 

    Key name:     scf_steps_offset
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     sqrt_interpolation
    Required:     no
    Key type:     boolean
    Expert:       Yes
    Experimental: No
    Default:      "false"
    Description:  Flag indicating whether or not to use square root technique for 
                  density interpolation. 

    Key name:     total_scf_steps_offset
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    0
    Max value:    0
    Default:      0
    Description:  

    Key name:     use_cmix
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "true"
    Description:  Use adaptive interpolation can improve energies/forces for some 
                  but not all systems. 

    Key name:     use_cpdgemr2d
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: No
    Default:      "true"
    Description:  if set to true, we use Cpdgemr2d to change matrix distribution 

    Key name:     use_gpu_fd
    Required:     no
    Key type:     boolean
    Expert:       No
    Experimental: Yes
    Default:      "false"
    Description:  Use gpus for kohn-sham orbital finite differencing. Depending on 
                  the balance of hardware characteristics this can provide a 
                  significant speedup but individual testing is required. 
                  Experimental. 

    Key name:     vxc_diag_nmax
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10000
    Default:      1
    Description:  Maximum band index for diagonal Vxc matrix elements. 

    Key name:     vxc_diag_nmin
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    1
    Max value:    10000
    Default:      1
    Description:  Minimum band index for diagonal Vxc matrix elements. 

    Key name:     wannier90_scdm
    Required:     no
    Key type:     integer
    Expert:       No
    Experimental: No
    Min value:    -2147483647
    Max value:    2
    Default:      0
    Description:  use scdm method to set the trial wannier functions 

    Key name:     wannier90_scdm_sigma
    Required:     no
    Key type:     double
    Expert:       No
    Experimental: No
    Min value:    0.000000e+00
    Max value:    unlimited
    Default:      1.000000
    Description:  when wannier90 is used to build wannier functions, the energy 
                  window parameter 

    Key name:     z_average_output_mode
    Required:     no
    Key type:     string
    Expert:       No
    Experimental: No
    Default:      "None"
    Allowed:      "wave functions" "potential and charge density" "None" 
    Description:  z_average_output_mode. 

 -- A Real Space Multigrid Electronic structure code --
 --      More information at www.rmgdft.org          --


[\[Introduction \] ](#introduction)
[\[Control options\] ](#control-options)
[\[Cell parameter options\] ](#cell-parameter-options)
[\[Pseudopotential related options\] ](#pseudopotential-related-options)
[\[Kohn Sham solver options\] ](#kohn-sham-solver-options)
[\[Exchange correlation options\] ](#exchange-correlation-options)
[\[Orbital occupation options\] ](#orbital-occupation-options)
[\[Charge density mixing options\] ](#charge-density-mixing-options)
[\[Relaxation and Molecular dynamics options\] ](#relaxation-and-molecular-dynamics-options)
[\[Diagonalization options\] ](#diagonalization-options)
[\[Performance related options\] ](#performance-related-options)
[\[LDAU options\] ](#ldau-options)
[\[TDDFT related options\] ](#tddft-related-options)
[\[Poisson solver options\] ](#poisson-solver-options)
[\[\] ](#)
[\[Testing options\] ](#testing-options)
[\[Miscellaneous options\] ](#miscellaneous-options)


## Introduction

    The RMG input file consists of a set of key-value pairs of the form.
    
        name = "scalar"

    
    where scalar can be an integer, double or boolean value.

        period_of_diagonalization = "1"
        charge_density_mixing = "0.5"
        initial_diagonalization = "true"
    
    There are also strings and arrays which are delineated by double quotes so an
    integer array with three elements would be.
    
        processor_grid = "2 2 2"
    
    while a string example would be
    
        description = "64 atom diamond cell test run at the gamma point"
    
    strings can span multiple lines so the following would be valid as well.
    
        description = "64 atom diamond cell test run at gamma point
        using a Vanderbilt ultrasoft pseudopotential"
    
    string vectors span multiple lines and are used to enter items like a kpoint list one item per line.
    
        kpoints = "0.00  0.00  0.00   0.50
                   0.25  0.25  0.50   0.25
                   0.50  0.50  0.50   0.25"


## Control options

<pre>
    <b>Key name:</b>     AFM
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  If true, anti-feromagnetic will be forced by symmetry operation if 
                  possible. 

    <b>Key name:</b>     BerryPhase
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  turn on/off Berry Phase calcualtion 

    <b>Key name:</b>     BerryPhaseCycle
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    10
    <b>Default:</b>      1
    <b>Description:</b>  Berry Phase loop without updating rho and potentials 

    <b>Key name:</b>     BerryPhaseDirection
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2
    <b>Default:</b>      2
    <b>Description:</b>  Berry Phase direction: it will be efield direction when efield is 
                  non zero 

    <b>Key name:</b>     STM_bias
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "-1.0 1.0"
    <b>Allowed:</b>      
    <b>Description:</b>  Bias (in unit of Volt) for STM calculation 

    <b>Key name:</b>     STM_height
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "2.0 4.0"
    <b>Allowed:</b>      
    <b>Description:</b>  Height range for STM calculation 

    <b>Key name:</b>     a_length
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  First lattice constant. 

    <b>Key name:</b>     adaptive_cmix
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    10.000000
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  Manual setting for the adaptive interpolation parameter. 

    <b>Key name:</b>     afd_cfac
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    10.000000
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  Manual setting for the adaptive finite differencing parameter. 

    <b>Key name:</b>     b_length
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  Second lattice constant. 

    <b>Key name:</b>     c_length
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  Third lattice constant. 

    <b>Key name:</b>     calculation_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Quench Electrons"
    <b>Allowed:</b>      "NSCF" "Constant Temperature And Energy" "Constant Volume And 
                  Energy" "Dimer Relax" "TDDFT" "Relax Structure" "Constant Pressure 
                  And Energy" "Quench Electrons" "Plot" "Psi Plot" "Band Structure 
                  Only" "STM" "NEB Relax" "Exx Only" 
    <b>Description:</b>  Type of calculation to perform. 

    <b>Key name:</b>     cell_relax
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  flag to control unit cell relaxation 

    <b>Key name:</b>     coalesce_factor
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    16
    <b>Default:</b>      4
    <b>Description:</b>  Grid coalescing factor. 

    <b>Key name:</b>     coalesce_states
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not to coalesce states. 

    <b>Key name:</b>     compressed_infile
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Flag indicating whether or not parallel restart wavefunction file 
                  uses compressed format. 

    <b>Key name:</b>     compressed_outfile
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Flag indicating whether or not parallel output wavefunction file 
                  uses compressed format. 

    <b>Key name:</b>     davidson_1stage_ortho
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Flag that improves davidson convergence for hard cases. Defaults 
                  to true but can be disabled for well behaved systems enabling 
                  higher performance. 

    <b>Key name:</b>     davidson_2stage_ortho
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag to futher improve davidson convergence for hard cases. 
                  Experimental. 

    <b>Key name:</b>     description
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  Description of the run. 

    <b>Key name:</b>     drho_precond_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Resta"
    <b>Allowed:</b>      "Kerker" "Resta" 
    <b>Description:</b>  Density mixing preconditioner method. Resta or Kerker are 
                  supported. 

    <b>Key name:</b>     energy_convergence_criterion
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000e-20
    <b>Max value:</b>    1.000000e-07
    <b>Default:</b>      1.000000e-10
    <b>Description:</b>  The RMS value of the estimated change in the total energy per step 
                  where we assume self consistency has been achieved. 

    <b>Key name:</b>     energy_output_units
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Hartrees"
    <b>Allowed:</b>      "Rydbergs" "Hartrees" 
    <b>Description:</b>  Units to be used when writing energy values to the output file. 
                  Hartrees or Rydbergs are available. 

    <b>Key name:</b>     epsg_guard
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000e-05
    <b>Default:</b>      1.000000e-07
    <b>Description:</b>  GGA guard value for low density regions. 

    <b>Key name:</b>     exx_integrals_filepath
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "afqmc_rmg"
    <b>Allowed:</b>      
    <b>Description:</b>  File/path for exact exchange integrals. 

    <b>Key name:</b>     exx_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Local fft"
    <b>Allowed:</b>      "Local fft" "Distributed fft" 
    <b>Description:</b>  FFT mode for exact exchange computations. 

    <b>Key name:</b>     exxdiv_treatment
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "gygi-baldereschi"
    <b>Allowed:</b>      "none" "gygi-baldereschi" 
    <b>Description:</b>  Exact exchange method for handling exx divergence at G=0. 

    <b>Key name:</b>     freeze_ldaU_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    2147483647
    <b>Default:</b>      500
    <b>Description:</b>  freeze the ldaU occupations ns_occ after this step. 

    <b>Key name:</b>     gpu_managed_memory
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Some AMD and Nvidia GPUs support managed gou memory which is 
                  useful when GPU memory limits are exceeded. 

    <b>Key name:</b>     input_tddft_file
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Waves/wave_tddft.out"
    <b>Allowed:</b>      
    <b>Description:</b>  Input file/path to read wavefunctions and other binary data from 
                  on a restart. 

    <b>Key name:</b>     input_wave_function_file
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Waves/wave.out"
    <b>Allowed:</b>      
    <b>Description:</b>  Input file/path to read wavefunctions and other binary data from 
                  on a restart. 

    <b>Key name:</b>     interpolation_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "FFT"
    <b>Allowed:</b>      "FFT" "prolong" "Cubic Polynomial" 
    <b>Description:</b>  Interpolation method for transferring data between the potential 
                  grid and the wavefunction grid. Mostly for diagnostic purposes. 

    <b>Key name:</b>     lambda_max
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000
    <b>Max value:</b>    100.000000
    <b>Default:</b>      3.000000
    <b>Description:</b>  Chebyshev smoothing parameter. Don't change unless you know what 
                  you're doing. 

    <b>Key name:</b>     lambda_min
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    2.000000
    <b>Default:</b>      0.300000
    <b>Description:</b>  Chebyshev smoothing parameter. Don't change unless you know what 
                  you're doing. 

    <b>Key name:</b>     max_exx_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    2147483647
    <b>Default:</b>      100
    <b>Description:</b>  Maximum number of self consistent steps to perform with hybrid 
                  functionals. 

    <b>Key name:</b>     max_scf_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      100
    <b>Description:</b>  Maximum number of self consistent steps to perform. Inner loop for 
                  hybrid functionals. 

    <b>Key name:</b>     noncollinear
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  If set true, noncollinear calculation. 

    <b>Key name:</b>     nvme_orbitals
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not orbitals should be mapped to disk. 

    <b>Key name:</b>     nvme_orbitals_filepath
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Orbitals/"
    <b>Allowed:</b>      
    <b>Description:</b>  File/path for runtime disk storage of orbitals. 

    <b>Key name:</b>     nvme_weights
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not projector weights should be mapped 
                  to disk. 

    <b>Key name:</b>     nvme_weights_filepath
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Weights/"
    <b>Allowed:</b>      
    <b>Description:</b>  File/path for disk storage of projector weights. 

    <b>Key name:</b>     nvme_work
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not work arrays should be mapped to 
                  disk. 

    <b>Key name:</b>     nvme_work_filepath
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Work/"
    <b>Allowed:</b>      
    <b>Description:</b>  File/path for disk storage of workspace. 

    <b>Key name:</b>     omp_threads_per_node
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    64
    <b>Default:</b>      0
    <b>Description:</b>  Number of Open MP threads each MPI process will use. A value of 0 
                  selects automatic setting. 

    <b>Key name:</b>     output_tddft_file
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Waves/wave_tddft.out"
    <b>Allowed:</b>      
    <b>Description:</b>  Output file/path to store wavefunctions and other binary data. 

    <b>Key name:</b>     output_wave_function_file
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Waves/wave.out"
    <b>Allowed:</b>      
    <b>Description:</b>  Output file/path to store wavefunctions and other binary data. 

    <b>Key name:</b>     pseudo_dir
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "."
    <b>Allowed:</b>      
    <b>Description:</b>  Directory where pseudopotentials are stored. 

    <b>Key name:</b>     qfunction_filepath
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Qfunctions/"
    <b>Allowed:</b>      
    <b>Description:</b>  File/path for runtime disk storage of qfunctions. 

    <b>Key name:</b>     qmc_nband
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      0
    <b>Description:</b>  The number of band used in rmg-qmcpack interface. 

    <b>Key name:</b>     read_serial_restart
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Directs RMG to read from serial restart files. Normally used when 
                  changing the sprocessor topology used during a restart run 

    <b>Key name:</b>     rms_convergence_criterion
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000e-03
    <b>Default:</b>      1.000000e-07
    <b>Description:</b>  The RMS value of the change in the total potential from step to 
                  step where we assume self consistency has been achieved. 

    <b>Key name:</b>     semilocal_projectors
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    6
    <b>Max value:</b>    40
    <b>Default:</b>      10
    <b>Description:</b>  Controls the number of semilocal projectors. 

    <b>Key name:</b>     spinorbit
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  If set true, spinorbit coupling calculation. 

    <b>Key name:</b>     start_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "LCAO Start"
    <b>Allowed:</b>      "Modified LCAO Start" "Restart TDDFT" "Start TDDFT" "FIREBALL 
                  Start" "Gaussian Start" "LCAO Start" "Restart From File" "Random 
                  Start" 
    <b>Description:</b>  Type of run. 

    <b>Key name:</b>     stress
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  flag to control stress cacluation 

    <b>Key name:</b>     stress_convergence_criterion
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    50.000000
    <b>Default:</b>      0.500000
    <b>Description:</b>  The stress criteria 

    <b>Key name:</b>     system_charge
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -unlimited
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  Number of excess holes in the system (useful for doped systems). 
                  Example: 2 means system is missing two electrons 

    <b>Key name:</b>     time_reversal
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  if false, no k -> -k symmetry 

    <b>Key name:</b>     use_energy_correction
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  Experimental energy correction term 

    <b>Key name:</b>     vdw_corr
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "None"
    <b>Allowed:</b>      "DFT-D3" "Grimme-D2" "DFT-D2" "None" 
    <b>Description:</b>  Type of vdw correction 

    <b>Key name:</b>     vdwdf_grid_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Fine"
    <b>Allowed:</b>      "Fine" "Coarse" 
    <b>Description:</b>  Type of grid to use when computing vdw-df correlation. 

    <b>Key name:</b>     vdwdf_kernel_filepath
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "vdW_kernel_table"
    <b>Allowed:</b>      
    <b>Description:</b>  File/path for vdW_kernel_table data. 

    <b>Key name:</b>     verbose
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag for writing out extra information 

    <b>Key name:</b>     wannier90
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  set up informations for wannier90 interface 

    <b>Key name:</b>     wannier90_scdm_mu
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -unlimited
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  when wannier90 is used to build wannier functions, the energy 
                  window parameter 

    <b>Key name:</b>     write_data_period
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -5
    <b>Max value:</b>    50000
    <b>Default:</b>      -1
    <b>Description:</b>  How often to write checkpoint files during the initial quench in 
                  units of SCF steps. During structural relaxations of molecular 
                  dynamics checkpoints are written each ionic step. 

    <b>Key name:</b>     write_qmcpack_restart
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  If true then a QMCPACK restart file is written as well as a serial 
                  restart file. 

    <b>Key name:</b>     write_qmcpack_restart_localized
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  If true then a QMCPACK restart file for localized orbitals 

    <b>Key name:</b>     write_serial_restart
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  RMG normally writes parallel restart files. These require that 
                  restarts have the same processor topology. If write_serial_restart 
                  = "true" then RMG will also write a serial restart file that can 
                  be used with a different processor topology 

</pre>
## Cell parameter options

<pre>
    <b>Key name:</b>     atomic_coordinate_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Absolute"
    <b>Allowed:</b>      "Absolute" "Cell Relative" 
    <b>Description:</b>  Flag indicated whether or not atomic coordinates are absolute or 
                  cell relative. 

    <b>Key name:</b>     bravais_lattice_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "None"
    <b>Allowed:</b>      "Orthorhombic Primitive" "Monoclinic Primitive" "Tetragonal 
                  Primitive" "Hexagonal Primitive" "Cubic Body Centered" "Triclinic 
                  Primitive" "Cubic Face Centered" "Cubic Primitive" "None" 
    <b>Description:</b>  Bravais Lattice Type. 

    <b>Key name:</b>     cell_movable
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "0 0 0 0 0 0 0 0 0 "
    <b>Description:</b>  9 numbers to control cell relaxation 

    <b>Key name:</b>     crds_units
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Bohr"
    <b>Allowed:</b>      "Angstrom" "Bohr" 
    <b>Description:</b>  Units for the atomic coordinates. 

    <b>Key name:</b>     frac_symmetry
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  For supercell calculation, one can disable the fractional 
                  translation symmetry 

    <b>Key name:</b>     grid_spacing
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.350000
    <b>Description:</b>  Approximate grid spacing (bohr). 

    <b>Key name:</b>     kpoint_distribution
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -2147483647
    <b>Max value:</b>    2147483647
    <b>Default:</b>      -1
    <b>Description:</b>  This option affects kpoint parallelization. If there are M MPI 
                  procs then N = M/kpoint_distribution procs are assigned to each 
                  kpoint. M must be evenly divisible by kpoint_distribution. 

    <b>Key name:</b>     kpoint_is_shift
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "0 0 0 "
    <b>Description:</b>  Three-D layout of the kpoint shift. 

    <b>Key name:</b>     kpoint_mesh
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "1 1 1 "
    <b>Description:</b>  Three-D layout of the kpoint mesh. 

    <b>Key name:</b>     kpoint_units
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Reciprocal lattice"
    <b>Allowed:</b>      "2pi/alat" "Reciprocal lattice" 
    <b>Description:</b>  kpoint units for reading kpoint 

    <b>Key name:</b>     lattice_units
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Bohr"
    <b>Allowed:</b>      "Alat" "Angstrom" "Bohr" 
    <b>Description:</b>  Units for the lattice vectors 

    <b>Key name:</b>     lattice_vector
    <b>Required:</b>     no
    <b>Key type:</b>     double array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "0 0 0 0 0 0 0 0 0 "
    <b>Description:</b>  The simulation cell may be specified using either lattice vectors, 
                  a0, a1, a2 or by lattice constants and a bravais lattice type. If 
                  lattice vectors are used they should be entered as a 3x3 matrix. 

    <b>Key name:</b>     ldos_end_grid
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "-1 -1 -1 "
    <b>Description:</b>  a ending grid point for ldos caclualtion 

    <b>Key name:</b>     ldos_start_grid
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "-1 -1 -1 "
    <b>Description:</b>  a grid point for starting ldos caclualtion 

    <b>Key name:</b>     potential_grid_refinement
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    4
    <b>Default:</b>      2
    <b>Description:</b>  Ratio of the potential grid density to the wavefunction grid 
                  density. For example if the wavefunction grid is (72,72,72) and 
                  potential_grid_refinement = "2" then the potential grid would be 
                  (144,144,144). The default value is 2 but it may sometimes be 
                  beneficial to adjust this. (For USPP the minimum value is also 2 
                  and it cannot be set lower. NCPP can be set to 1). 

    <b>Key name:</b>     processor_grid
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "1 1 1 "
    <b>Description:</b>  Three-D (x,y,z) layout of the MPI processes. 

    <b>Key name:</b>     sts_end_grid
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "-1 -1 -1 "
    <b>Description:</b>  a ending grid point for sts caclualtion 

    <b>Key name:</b>     sts_start_grid
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "-1 -1 -1 "
    <b>Description:</b>  a grid point for starting sts caclualtion 

    <b>Key name:</b>     use_symmetry
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2
    <b>Default:</b>      2
    <b>Description:</b>  0: never use symmetry, 1: always use symmetry, 

    <b>Key name:</b>     wavefunction_grid
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "1 1 1 "
    <b>Description:</b>  Three-D (x,y,z) dimensions of the grid the wavefunctions are 
                  defined on. 

</pre>
## Pseudopotential related options

<pre>
    <b>Key name:</b>     all_electron_parm
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> Yes
    <b>Min value:</b>    1
    <b>Max value:</b>    12
    <b>Default:</b>      4
    <b>Description:</b>  Gygi all electron parameter. 

    <b>Key name:</b>     atomic_orbital_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "delocalized"
    <b>Allowed:</b>      "delocalized" "localized" 
    <b>Description:</b>  Atomic Orbital Type. Choices are localized and delocalized. 

    <b>Key name:</b>     energy_cutoff_parameter
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.600000
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.800000
    <b>Description:</b>  

    <b>Key name:</b>     filter_factor
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.060000
    <b>Max value:</b>    1.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  Filtering factor. 

    <b>Key name:</b>     internal_pseudo_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "sg15"
    <b>Allowed:</b>      "all_electron" "nc_standard" "nc_accuracy" "sg15" "ultrasoft" 
    <b>Description:</b>  Internal pseudopotential type. Choices are sg15, ultrasoft, 
                  nc_accuracy or all_electron 

    <b>Key name:</b>     localize_localpp
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  The local potential associated with a particular ion also decays 
                  rapidly in real-space with increasing r. As with beta projectors 
                  truncating the real-space representation for large cells can lead 
                  to significant computational savings with a small loss of accuracy 
                  but it should be set to false for small cells. 

    <b>Key name:</b>     localize_projectors
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  The Beta function projectors for a particular ion decay rapidly in 
                  real-space with increasing r. For large cells truncating the 
                  real-space representation of the projector can lead to significant 
                  computational savings with a small loss of accuracy. For smaller 
                  cells the computational cost is the same for localized or 
                  delocalized projectors so it is better to set localize_projectors 
                  to false. 

    <b>Key name:</b>     max_nlradius
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    2.000000
    <b>Max value:</b>    10000.000000
    <b>Default:</b>      10000.000000
    <b>Description:</b>  maximum radius for non-local projectors 

    <b>Key name:</b>     max_qradius
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    2.000000
    <b>Max value:</b>    10000.000000
    <b>Default:</b>      10000.000000
    <b>Description:</b>  maximum radius for qfunc in ultra-pseudopotential 

    <b>Key name:</b>     min_nlradius
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000
    <b>Max value:</b>    10000.000000
    <b>Default:</b>      2.000000
    <b>Description:</b>  minimum radius for non-local projectors 

    <b>Key name:</b>     min_qradius
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000
    <b>Max value:</b>    10000.000000
    <b>Default:</b>      2.000000
    <b>Description:</b>  minimum radius for qfunc in ultra-pseudopotential 

    <b>Key name:</b>     projector_expansion_factor
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.500000
    <b>Max value:</b>    3.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  When using localized projectors the radius can be adjusted with 
                  this parameter. 

    <b>Key name:</b>     pseudopotential
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  External pseudopotentials may be specfied with this input key. The 
                  format uses the atomic symbol followed by the pseudopotential file 
                  name. pseudopotential = "Ni Ni.UPF O O.UPF" 

    <b>Key name:</b>     use_bessel_projectors
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  When a semi-local pseudopotential is being used projectors will be 
                  generated using Bloechl's procedure with Bessel functions as the 
                  basis set if this is true. 

</pre>
## Kohn Sham solver options

<pre>
    <b>Key name:</b>     davidson_max_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    5
    <b>Max value:</b>    20
    <b>Default:</b>      8
    <b>Description:</b>  Maximum number of iterations for davidson diagonalization. 

    <b>Key name:</b>     davidson_multiplier
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    6
    <b>Default:</b>      0
    <b>Description:</b>  The davidson solver expands the eigenspace with the maximum 
                  expansion factor being set by the value of davidson_multiplier. 
                  Larger values often lead to faster convergence but because the 
                  computational cost of the davidson diagonalization step scales as 
                  the cube of the number of eigenvectors the optimal value based on 
                  the fastest time to solution depends on the number of orbitals. If 
                  not specified explicitly or set to 0 RMG uses the following 
                  algorithm to set the value. 
                  
                  Number of orbitals <= 600 davidson_multiplier= "4" 
                  600 < Number of orbitals <= 900 davidson_multiplier = "3" 
                  Number of orbitals > 900 davidson_multiplier = "2" 
                  
                  For very large problems the N^3 scaling makes even a factor of 2 
                  prohibitively costly and the multigrid solver is a better choice. 

    <b>Key name:</b>     davidson_premg
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    8
    <b>Default:</b>      0
    <b>Description:</b>  If the davidson solver is selected this parameter controls the 
                  number of multigrid steps to use before enabling davidson. 

    <b>Key name:</b>     kohn_sham_coarse_time_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.200000
    <b>Default:</b>      1.000000
    <b>Description:</b>  Time step to use in the kohn-sham multigrid solver on the coarse 
                  levels. 

    <b>Key name:</b>     kohn_sham_fd_order
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    6
    <b>Max value:</b>    12
    <b>Default:</b>      8
    <b>Description:</b>  RMG uses finite differencing to represent the kinetic energy 
                  operator and the accuracy of the representation is controllable by 
                  the kohn_sham_fd_order parameter. The default is 8 and is fine for 
                  most purposes but higher accuracy is obtainable with 10th or 12th 
                  order at the cost of some additional computational expense. 

    <b>Key name:</b>     kohn_sham_mg_levels
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -1
    <b>Max value:</b>    6
    <b>Default:</b>      -1
    <b>Description:</b>  Number of multigrid levels to use in the kohn-sham multigrid 
                  preconditioner. 

    <b>Key name:</b>     kohn_sham_mg_timestep
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    2.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  timestep for multigrid correction. 

    <b>Key name:</b>     kohn_sham_mucycles
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    6
    <b>Default:</b>      3
    <b>Description:</b>  Number of mu (also known as W) cycles to use in the kohn-sham 
                  multigrid preconditioner. 

    <b>Key name:</b>     kohn_sham_post_smoothing
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    5
    <b>Default:</b>      1
    <b>Description:</b>  Number of global grid post-smoothing steps to perform after a 
                  multigrid preconditioner iteration. 

    <b>Key name:</b>     kohn_sham_pre_smoothing
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    5
    <b>Default:</b>      2
    <b>Description:</b>  Number of global grid pre-smoothing steps to perform before a 
                  multigrid preconditioner iteration. 

    <b>Key name:</b>     kohn_sham_solver
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "davidson"
    <b>Allowed:</b>      "davidson" "multigrid" 
    <b>Description:</b>  RMG supports a pure multigrid Kohn-Sham solver as well as a 
                  multigrid preconditioned davidson solver. The davidson solver is 
                  usually better for smaller problems with the pure multigrid solver 
                  often being a better choice for very large problems. 

    <b>Key name:</b>     kohn_sham_time_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    2.000000
    <b>Default:</b>      0.660000
    <b>Description:</b>  Smoothing timestep to use on the fine grid in the the kohn-sham 
                  multigrid preconditioner. 

    <b>Key name:</b>     prolong_order
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    12
    <b>Default:</b>      10
    <b>Description:</b>  Debug option that controls interpolation order used to form the 
                  charge density and to compute the kinetic component of stress. If 
                  a value of 0 is selected then an FFT will be used. 

    <b>Key name:</b>     unoccupied_tol_factor
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000
    <b>Max value:</b>    100000.000000
    <b>Default:</b>      1000.000000
    <b>Description:</b>  When using the Davidson Kohn-Sham solver unoccupied states are 
                  converged to a less stringent tolerance than occupied orbitals 
                  with the ratio set by this parameter. 

    <b>Key name:</b>     use_block_diag
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not to use block diagonalization. 

    <b>Key name:</b>     use_rmm_diis
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Flag indicating whether or not to use the RMM-DIIS algorithm in 
                  the mulgrid solver. 

</pre>
## Exchange correlation options

<pre>
    <b>Key name:</b>     exchange_correlation_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "AUTO_XC"
    <b>Allowed:</b>      "hartree-fock" "VDW-DF" "vdw-df" "gaupbe" "HSE" "PBE0" "mgga tb09" 
                  "TB09" "VDW-DF-CX" "optbk88" "AUTO_XC" "BLYP" "hse" "tb09" "GGA 
                  BLYP" "GGA XB CP" "ev93" "GGA XP CP" "Q2D" "pbe0" "PW91" "pbe" 
                  "bp" "vdw-df-c09" "optb86b" "pbesol" "pw91" "b86bpbe" "BP" "blyp" 
                  "PBE" "pz" "PZ" "LDA" "hcth" "revpbe" "m06l" "REVPBE" "vdw-df-cx" 
                  "pw86pbe" "B3LYP" "PW86PBE" "PBESOL" "sla+pw+pbe+vdw1" "sogga" 
                  "q2d" "MGGA TB09" "GGA PBE" "wc" "tpss" "HCTH" "b3lyp" "olyp" 
    <b>Description:</b>  Most pseudopotentials specify the exchange correlation type they 
                  were generated with and the default value of AUTO_XC means that 
                  the type specified in the pseudopotial is what RMG will use. That 
                  can be overridden by specifying a value here. 

    <b>Key name:</b>     exx_convergence_criterion
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000e-12
    <b>Max value:</b>    1.000000e-06
    <b>Default:</b>      1.000000e-09
    <b>Description:</b>  Convergence criterion for the EXX delta from step to step where we 
                  assume EXX consistency has been achieved. 

    <b>Key name:</b>     exx_fraction
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -1.000000
    <b>Max value:</b>    1.000000
    <b>Default:</b>      -1.000000
    <b>Description:</b>  when hybrid functional is used, the fraction of Exx 

    <b>Key name:</b>     vexx_fft_threshold
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000e-14
    <b>Max value:</b>    0.100000
    <b>Default:</b>      1.000000e-14
    <b>Description:</b>  The value for the EXX delta where we switch from single to double 
                  precision ffts. Single precision ffts are generally accurate 
                  enough. 

    <b>Key name:</b>     x_gamma_extrapolation
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  if set true, use exx extrapolation to gamma 

</pre>
## Orbital occupation options

<pre>
    <b>Key name:</b>     MP_order
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    5
    <b>Default:</b>      2
    <b>Description:</b>  Order of Methefessel Paxton occupation. 

    <b>Key name:</b>     dos_broading
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.100000
    <b>Description:</b>  For DOS with Gaussian broading method 

    <b>Key name:</b>     dos_method
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "tetrahedra"
    <b>Allowed:</b>      "Gaussian" "tetrahedra" 
    <b>Description:</b>  tetrahedra or gauss smearing method for DOS calculation 

    <b>Key name:</b>     occupation_electron_temperature_eV
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    2.000000
    <b>Default:</b>      0.040000
    <b>Description:</b>  Target electron temperature when not using fixed occupations. 

    <b>Key name:</b>     occupation_number_mixing
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  Mixing parameter for orbital occupations when not using fixed 
                  occupations. 

    <b>Key name:</b>     occupations_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Fermi Dirac"
    <b>Allowed:</b>      "Tetrahedron" "MethfesselPaxton" "Cold Smearing" "Error Function" 
                  "Gaussian" "Fermi Dirac" "Fixed" 
    <b>Description:</b>  RMG supports several different ways of specifying orbital 
                  occupations. For a spin polarized system one may specify the 
                  occupations for up and down separately. In the case of a non-zero 
                  electronic temperature these will be adjusted as the calculation 
                  proceeds based on this setting. 

    <b>Key name:</b>     states_count_and_occupation
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  Occupation string for states. Format for a system with 240 
                  electrons and 20 unoccupied states would be. "120 2.0 20 0.0" 

    <b>Key name:</b>     states_count_and_occupation_spin_down
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  Occupation string for spin down states. Format is the same as for 
                  states_count_and_occupation. Total number of states must match 
                  spin up occupation string. 

    <b>Key name:</b>     states_count_and_occupation_spin_up
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  Occupation string for spin up states. Format is the same as for 
                  states_count_and_occupation. Total number of states must match 
                  spin down occupation string. 

    <b>Key name:</b>     tetra_method
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Bloechl"
    <b>Allowed:</b>      "Optimized" "Linear" "Bloechl" 
    <b>Description:</b>  tetrahedron method to use 

    <b>Key name:</b>     unoccupied_states_per_kpoint
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -2147483647
    <b>Max value:</b>    2147483647
    <b>Default:</b>      -1
    <b>Description:</b>  The number of unoccupied orbitals. A value that is 15-20% of the 
                  number of occupied orbitals generally works well. 

</pre>
## Charge density mixing options

<pre>
    <b>Key name:</b>     charge_broyden_order
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    12
    <b>Default:</b>      10
    <b>Description:</b>  Number of previous steps to use when Broyden mixing is used to 
                  update the charge density. 

    <b>Key name:</b>     charge_density_mixing
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.500000
    <b>Description:</b>  Proportion of the current charge density to replace with the new 
                  density after each scf step when linear mixing is used. 

    <b>Key name:</b>     charge_mixing_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Auto"
    <b>Allowed:</b>      "Broyden" "Auto" "Pulay" "Linear" 
    <b>Description:</b>  RMG supports Broyden, Pulay and Linear mixing When the davidson 
                  Kohn-Sham solver is selected Broyden or Pulay are preferred. For 
                  the multigrid solver Broyden is usually the best choice. 

    <b>Key name:</b>     charge_pulay_Gspace
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  if set true, charge density mixing the residual in G space 

    <b>Key name:</b>     charge_pulay_order
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    10
    <b>Default:</b>      5
    <b>Description:</b>  Number of previous steps to use when Pulay mixing is used to 
                  update the charge density. 

    <b>Key name:</b>     charge_pulay_refresh
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    2147483647
    <b>Default:</b>      100
    <b>Description:</b>  charge Pulay mixing reset steps. 

    <b>Key name:</b>     drho_precond
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  if set true, charge density residual is preconded with 
                  q^2/(q^2+q0^2) 

    <b>Key name:</b>     drho_precond_q0
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    10.000000
    <b>Default:</b>      0.250000
    <b>Description:</b>  Kerker type preconditioning the charge density residual by 
                  q^2/(q^2+q0^2) See Kresse and Furthmueller, Computational 
                  Materials Science 6 (1996) 15-50 

    <b>Key name:</b>     ldau_mixing
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  Proportion of the current ldau occupation to replace with the new 
                  ones after each scf step when linear mixing is used. 

    <b>Key name:</b>     ldau_mixing_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Linear"
    <b>Allowed:</b>      "Broyden" "Auto" "Pulay" "Linear" 
    <b>Description:</b>  RMG supports Pulay and Linear mixing for DFT+U occupation mixing 

    <b>Key name:</b>     ldau_pulay_order
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    10
    <b>Default:</b>      5
    <b>Description:</b>  Number of previous steps to use when Pulay mixing is used to 
                  update the ldau occupation . 

    <b>Key name:</b>     ldau_pulay_refresh
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    2147483647
    <b>Default:</b>      100
    <b>Description:</b>  ldau pulay mixing reset steps 

    <b>Key name:</b>     ldau_pulay_scale
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  

    <b>Key name:</b>     potential_acceleration_constant_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    4.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  When set to a non-zero value this parameter causes RMG to perform 
                  a band by band update of the self-consistent potential during the 
                  course of an SCF step when the multigrid kohn_sham_solver is 
                  chosen. This means that updates to the lower energy orbitals are 
                  incorporated into the SCF potential seen by the higher energy 
                  orbitals as soon as they are computed. This can lead to faster 
                  convergence and better stability for many systems. The option 
                  should only be used with Linear mixing. Even when the davidson 
                  solver is chosen this parameter may be used since the first few 
                  steps with davidson usually uses the multigrid solver. 

    <b>Key name:</b>     resta_beta
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1000.000000
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  Beta parameter for resta charge density preconditioning. The 
                  default value of 0.0 means the value should be autmatically 
                  determined. 

</pre>
## Relaxation and Molecular dynamics options

<pre>
    <b>Key name:</b>     dynamic_time_counter
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    0
    <b>Default:</b>      0
    <b>Description:</b>  

    <b>Key name:</b>     dynamic_time_delay
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    5
    <b>Max value:</b>    5
    <b>Default:</b>      5
    <b>Description:</b>  

    <b>Key name:</b>     force_grad_order
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    12
    <b>Default:</b>      8
    <b>Description:</b>  Atomic forces may be computed to varying degrees of accuracy 
                  depending on the requirements of a specific problem. A value of 0 
                  implies highest accuracy which is obtained by using FFTs in place 
                  of finite differencing. 

    <b>Key name:</b>     ionic_time_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      50.000000
    <b>Description:</b>  Ionic time step for use in molecular dynamics and structure 
                  optimizations. 

    <b>Key name:</b>     ionic_time_step_decrease
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.500000
    <b>Description:</b>  Factor by which ionic timestep is decreased when dynamic timesteps 
                  are enabled. 

    <b>Key name:</b>     ionic_time_step_increase
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000
    <b>Max value:</b>    3.000000
    <b>Default:</b>      1.100000
    <b>Description:</b>  Factor by which ionic timestep is increased when dynamic timesteps 
                  are enabled. 

    <b>Key name:</b>     max_ionic_time_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    150.000000
    <b>Default:</b>      150.000000
    <b>Description:</b>  Maximum ionic time step to use for molecular dynamics or 
                  structural optimizations. 

    <b>Key name:</b>     max_md_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      100
    <b>Description:</b>  Maximum number of molecular dynamics steps to perform. 

    <b>Key name:</b>     md_integration_order
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "5th Beeman-Velocity Verlet"
    <b>Allowed:</b>      "5th Beeman-Velocity Verlet" "3rd Beeman-Velocity Verlet" "2nd 
                  Velocity Verlet" 
    <b>Description:</b>  Integration order for molecular dynamics. 

    <b>Key name:</b>     md_nose_oscillation_frequency_THz
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      15.590000
    <b>Description:</b>  

    <b>Key name:</b>     md_number_of_nose_thermostats
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    5
    <b>Max value:</b>    5
    <b>Default:</b>      5
    <b>Description:</b>  Number of Nose thermostats to use during Constant Volume and 
                  Temperature MD. 

    <b>Key name:</b>     md_randomize_velocity
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  The initial ionic velocities for a molecular dyanamics run are 
                  randomly initialized to the target temperature. 

    <b>Key name:</b>     md_temperature
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      300.000000
    <b>Description:</b>  Target MD Temperature. 

    <b>Key name:</b>     md_temperature_control
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Nose Hoover Chains"
    <b>Allowed:</b>      "Anderson Rescaling" "Nose Hoover Chains" 
    <b>Description:</b>  Type of temperature control method to use in molecular dynamics. 

    <b>Key name:</b>     relax_dynamic_timestep
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not to use dynamic timesteps in 
                  relaxation mode. 

    <b>Key name:</b>     relax_mass
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Atomic"
    <b>Allowed:</b>      "Equal" "Atomic" 
    <b>Description:</b>  Mass to use for structural relaxation, either atomic masses, or 
                  the mass of carbon for all atoms. 

    <b>Key name:</b>     relax_max_force
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      2.500000e-03
    <b>Description:</b>  Force value at which an ionic relaxation is considered to be 
                  converged. 

    <b>Key name:</b>     relax_method
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "LBFGS"
    <b>Allowed:</b>      "LBFGS" "MD Min" "Quick Min" "FIRE" "Fast Relax" 
    <b>Description:</b>  Type of relaxation method to use for structural optimizations. 

    <b>Key name:</b>     renormalize_forces
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Flag indicating whether or not to renormalize forces. 

    <b>Key name:</b>     tddft_frequency
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.200000
    <b>Description:</b>  TDDFT frequency for use in TDDFT vector potential mode 

    <b>Key name:</b>     tddft_time_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.200000
    <b>Description:</b>  TDDFT time step for use in TDDFT mode 

</pre>
## Diagonalization options

<pre>
    <b>Key name:</b>     extra_random_lcao_states
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      0
    <b>Description:</b>  LCAO (Linear Combination of Atomic Orbitals) is the default 
                  startup method for RMG. The atomic orbitals are obtained from the 
                  pseudpotentials but in some cases better convergence may be 
                  obtained by adding extra random wavefunctions in addition to the 
                  atomic orbitals. 

    <b>Key name:</b>     folded_spectrum
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  When the number of eigenvectors is large using folded_spectrum is 
                  substantially faster than standard diagonalization. It also tends 
                  to converge better for metallic systems. It works with the 
                  multigrid kohn_sham_solver but not the davidson solver. 

    <b>Key name:</b>     folded_spectrum_iterations
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    20
    <b>Default:</b>      2
    <b>Description:</b>  Number of folded spectrum iterations to perform. 

    <b>Key name:</b>     folded_spectrum_width
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.100000
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.300000
    <b>Description:</b>  Submatrix width to use as a fraction of the full spectrum. The 
                  folded spectrum width ranges from 0.10 to 1.0. For insulators and 
                  semiconductors a value of 0.3 is appropriate. For metals values 
                  between 0.15 to 0.2 tend to be better. The default value is 0.3 

    <b>Key name:</b>     initial_diagonalization
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Perform initial subspace diagonalization. 

    <b>Key name:</b>     period_of_diagonalization
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      1
    <b>Description:</b>  Diagonalization period (per scf step). Mainly for debugging and 
                  should not be changed for production. 

    <b>Key name:</b>     scalapack_block_factor
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    4
    <b>Max value:</b>    2147483647
    <b>Default:</b>      64
    <b>Description:</b>  Block size to use with scalapack. Optimal value is dependent on 
                  matrix size and system hardware. 

    <b>Key name:</b>     subdiag_driver
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "auto"
    <b>Allowed:</b>      "auto" "rocsolver" "elpa" "cusolver" "magma" "scalapack" "lapack" 
    <b>Description:</b>  Driver type used for subspace diagonalization of the eigenvectors. 

    <b>Key name:</b>     subdiag_groups
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    16
    <b>Default:</b>      1
    <b>Description:</b>  Number of scalapack or elpa groups. 

</pre>
## Performance related options

<pre>
    <b>Key name:</b>     fd_allocation_limit
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1024
    <b>Max value:</b>    262144
    <b>Default:</b>      65536
    <b>Description:</b>  Allocation sizes in finite difference routines less than this 
                  value are stack rather than heap based. 

    <b>Key name:</b>     mpi_queue_mode
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Use mpi queue mode. 

    <b>Key name:</b>     non_local_block_size
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    64
    <b>Max value:</b>    40000
    <b>Default:</b>      512
    <b>Description:</b>  Block size to use when applying the non-local and S operators. A 
                  value at least as large as the number of wavefunctions produces 
                  better performance but requires more memory. 

    <b>Key name:</b>     preconditioner_threshold
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000e-09
    <b>Max value:</b>    0.100000
    <b>Default:</b>      0.100000
    <b>Description:</b>  The RMS value of the change in the total potential where we switch 
                  the preconditioner from single to double precision. 

    <b>Key name:</b>     require_huge_pages
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  If set RMG assumes that sufficient huge pages are available. Bad 
                  things may happen if this is not true. 

    <b>Key name:</b>     rmg_threads_per_node
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    64
    <b>Default:</b>      0
    <b>Description:</b>  Number of Multigrid/Davidson threads each MPI process will use. A 
                  value of 0 means set automatically. 

    <b>Key name:</b>     spin_manager_thread
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  When mpi_queue_mode is enabled the manager thread spins instead of 
                  sleeping. 

    <b>Key name:</b>     spin_worker_threads
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  When mpi_queue_mode is enabled the worker threads spin instead of 
                  sleeping. 

    <b>Key name:</b>     state_block_size
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    2147483647
    <b>Default:</b>      64
    <b>Description:</b>  State block size used in nlforce. Larger values require more 
                  memory but can 

    <b>Key name:</b>     use_alt_zgemm
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not to use alternate zgemm 
                  implementation. 

    <b>Key name:</b>     use_async_allreduce
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  RMG uses MPI_Allreduce function calls in several places and for 
                  large problems these can account for a significant fraction of the 
                  total run time. In most cases using the asynchronous MPI versions 
                  of the functions is faster but this is not true for all platforms 
                  and in that casesetting this flag to false can improve 
                  performance. 

    <b>Key name:</b>     use_hwloc
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Use internal hwloc setup if available. If both this and use_numa 
                  are true hwloc takes precedence. 

    <b>Key name:</b>     use_numa
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Numa stands for Non Uniform Memory Access and means that the main 
                  memory of a computer is organized into seperate distinct banks. 
                  Each bank is then attached to a CPU core or group of cores and 
                  while all cores can normally access all banks the access speed is 
                  faster for directly attached banks. Ensuring that individual CPU 
                  cores mostly access memory in banks they are directly attached to 
                  can have a large impact on performance. Process mapping that does 
                  this can normally be done when jobs are submitted and run via 
                  arguments to mpirun/mpiexec but if this is not done RMG will 
                  attempt to provide an optimal mapping if use_numa is set to true. 

</pre>
## LDAU options

<pre>
    <b>Key name:</b>     Hubbard_U
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  Hubbard U parameter for each atomic species using the format 
                  Hubbard_U="Ni 6.5 3d 0.0 0.0 0.0" 

    <b>Key name:</b>     ldaU_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "None"
    <b>Allowed:</b>      "Simple" "None" 
    <b>Description:</b>  Type of lda+u implementation. 

    <b>Key name:</b>     ldaU_radius
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000
    <b>Max value:</b>    12.000000
    <b>Default:</b>      9.000000
    <b>Description:</b>  Max radius of atomic orbitals to be used in LDA+U projectors. 

</pre>
## TDDFT related options

<pre>
    <b>Key name:</b>     restart_tddft
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  restart TDDFT 

    <b>Key name:</b>     tddft_gpu
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  use gpu for ELYDYN or not 

    <b>Key name:</b>     tddft_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "electric field"
    <b>Allowed:</b>      "vector potential" "point charge" "electric field" 
    <b>Description:</b>  TDDFT mode 

    <b>Key name:</b>     tddft_noscf
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  TDDFT run read data directly from the last scf job 

    <b>Key name:</b>     tddft_qgau
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      1.000000
    <b>Description:</b>  Gaussian parameter for point charge to Gaussian charge 

    <b>Key name:</b>     tddft_qpos
    <b>Required:</b>     no
    <b>Key type:</b>     double array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "0 0 0 "
    <b>Description:</b>  cartesian coordinate of the point charge for tddft 

    <b>Key name:</b>     tddft_start_state
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      0
    <b>Description:</b>  the starting state to use in tddft dynamics 

    <b>Key name:</b>     tddft_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      2000
    <b>Description:</b>  Maximum number of tddft steps to perform. 

</pre>
## Poisson solver options

<pre>
    <b>Key name:</b>     hartree_max_sweeps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    5
    <b>Max value:</b>    100
    <b>Default:</b>      10
    <b>Description:</b>  Maximum number of hartree iterations to perform per scf step. 

    <b>Key name:</b>     hartree_min_sweeps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    5
    <b>Default:</b>      5
    <b>Description:</b>  Minimum number of hartree iterations to perform per scf step. 

    <b>Key name:</b>     hartree_rms_ratio
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1000.000000
    <b>Max value:</b>    unlimited
    <b>Default:</b>      100000.000000
    <b>Description:</b>  Ratio between target RMS for get_vh and RMS total potential. 

    <b>Key name:</b>     poisson_coarse_time_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.400000
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.800000
    <b>Description:</b>  Time step to use in the poisson multigrid solver on the coarse 
                  levels. 

    <b>Key name:</b>     poisson_coarsest_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    10
    <b>Max value:</b>    100
    <b>Default:</b>      25
    <b>Description:</b>  Number of smoothing steps to use on the coarsest level in the 
                  hartree multigrid solver. 

    <b>Key name:</b>     poisson_finest_time_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.400000
    <b>Max value:</b>    1.000000
    <b>Default:</b>      1.000000
    <b>Description:</b>  Time step to use in the poisson multigrid solver on the finest 
                  level. 

    <b>Key name:</b>     poisson_mg_levels
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -1
    <b>Max value:</b>    6
    <b>Default:</b>      -1
    <b>Description:</b>  Number of multigrid levels to use in the hartree multigrid solver. 

    <b>Key name:</b>     poisson_mucycles
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    4
    <b>Default:</b>      3
    <b>Description:</b>  Number of mu (also known as W) cycles to use in the hartree 
                  multigrid solver. 

    <b>Key name:</b>     poisson_post_smoothing
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    6
    <b>Default:</b>      1
    <b>Description:</b>  Number of global hartree grid post-smoothing steps to perform 
                  after a multigrid iteration. 

    <b>Key name:</b>     poisson_pre_smoothing
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    6
    <b>Default:</b>      2
    <b>Description:</b>  Number of global hartree grid pre-smoothing steps to perform 
                  before a multigrid iteration. 

    <b>Key name:</b>     poisson_solver
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "pfft"
    <b>Allowed:</b>      "pfft" "multigrid" 
    <b>Description:</b>  poisson solver. 

</pre>
## 

<pre>
    <b>Key name:</b>     charge_analysis_period
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    500
    <b>Default:</b>      10
    <b>Description:</b>  How often to perform and write out charge analysis. 

    <b>Key name:</b>     output_rho_xsf
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Generate xsf format for electronic density. 

    <b>Key name:</b>     write_eigvals_period
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    100
    <b>Default:</b>      5
    <b>Description:</b>  How often to output eigenvalues in units of scf steps. 

    <b>Key name:</b>     write_orbital_overlaps
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  If true the orbital overlap matrix from successive MD steps is 
                  written. 

    <b>Key name:</b>     write_pdos
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag to write partial density of states. 

    <b>Key name:</b>     write_pseudopotential_plots
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag to indicate whether or not to write pseudopotential plots. 

</pre>
## Testing options

<pre>
    <b>Key name:</b>     test_bond_length
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    20.000000
    <b>Default:</b>      nan
    <b>Description:</b>  Expected dimer bond length for testing. 

    <b>Key name:</b>     test_bond_length_tolerance
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000e-04
    <b>Max value:</b>    0.100000
    <b>Default:</b>      1.000000e-03
    <b>Description:</b>  Test bond length tolerance. 

    <b>Key name:</b>     test_energy
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -1000000000.000000
    <b>Max value:</b>    1000000000.000000
    <b>Default:</b>      nan
    <b>Description:</b>  Expected final energy for testing. 

    <b>Key name:</b>     test_energy_tolerance
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1.000000e-08
    <b>Max value:</b>    1.000000e-04
    <b>Default:</b>      1.000000e-07
    <b>Description:</b>  Test final energy tolerance. 

    <b>Key name:</b>     test_steps
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    1000
    <b>Default:</b>      0
    <b>Description:</b>  Expected number of scf steps for testing. 

    <b>Key name:</b>     test_steps_tolerance
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    10
    <b>Default:</b>      1
    <b>Description:</b>  Test scf steps tolerance. 

</pre>
## Miscellaneous options

<pre>
    <b>Key name:</b>     E_POINTS
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    2147483647
    <b>Default:</b>      201
    <b>Description:</b>  number of E points for ldos and sts calculation 

    <b>Key name:</b>     Emax
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -100.000000
    <b>Max value:</b>    100.000000
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  

    <b>Key name:</b>     Emin
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -100.000000
    <b>Max value:</b>    100.000000
    <b>Default:</b>      -6.000000
    <b>Description:</b>  

    <b>Key name:</b>     ExxCholMax
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    64
    <b>Default:</b>      8
    <b>Description:</b>  maximum number of Exx integral cholesky vectors 

    <b>Key name:</b>     ExxIntCholosky
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  if set true, Exx integrals are Cholesky factorized to 3-index 

    <b>Key name:</b>     adaptive_convergence
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Parameters that control initial SCF convergence are adaptively 
                  modified. Affected parameters include density mixing and 
                  preconditioning and electron temperature. 

    <b>Key name:</b>     alt_laplacian
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Flag indicating whether or not to use alternate laplacian weights 
                  for some operators. 

    <b>Key name:</b>     boundary_condition_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Periodic"
    <b>Allowed:</b>      "Periodic" 
    <b>Description:</b>  Boundary condition type Only periodic is currently implemented. 

    <b>Key name:</b>     charge_analysis
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Voronoi"
    <b>Allowed:</b>      "Voronoi" "None" 
    <b>Description:</b>  Type of charge analysis to use. Only Voronoi deformation density 
                  is currently available. 

    <b>Key name:</b>     cube_pot
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  if set true, total potential is printed out in cube format 

    <b>Key name:</b>     cube_rho
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  if set true, charge density is printed out in cube format 

    <b>Key name:</b>     cube_states_list
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  plot the states listed here 

    <b>Key name:</b>     cube_vh
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  if set true, hatree potential is printed out in cube format 

    <b>Key name:</b>     dftd3_version
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    2
    <b>Max value:</b>    6
    <b>Default:</b>      3
    <b>Description:</b>  Grimme's DFT-D3 versions, 

    <b>Key name:</b>     dipole_correction
    <b>Required:</b>     no
    <b>Key type:</b>     integer array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "0 0 0 "
    <b>Description:</b>  (1,1,1) for molecule, dipole correction in all directions. 

    <b>Key name:</b>     dipole_moment
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Turns on calculation of dipole moment for the entire cell. 

    <b>Key name:</b>     ecutrho
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    10000.000000
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  ecut for rho in unit of Ry. 

    <b>Key name:</b>     ecutwfc
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    10000.000000
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  ecut for wavefunctions in unit of Ry. 

    <b>Key name:</b>     electric_field
    <b>Required:</b>     no
    <b>Key type:</b>     double array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "0 0 0 "
    <b>Description:</b>  Components of the electric field in reciprocal lattice direction 
                  and unit of Ha/bohr. 

    <b>Key name:</b>     electric_field_tddft
    <b>Required:</b>     no
    <b>Key type:</b>     double array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "0 0 0 "
    <b>Description:</b>  the electric field for TDDFT in reciprocal lattice direction and 
                  unit of Ha/bohr 

    <b>Key name:</b>     equal_initial_density
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Specifies whether to set initial up and down density to be equal. 

    <b>Key name:</b>     exx_int_flag
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  If set true, calculate the exact exchange integrals. 

    <b>Key name:</b>     fast_density
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Use a faster but less accurate method to generate the charge 
                  density from the electronic wavefunctions. As the cutoff 
                  (grid-density) increases this method improves in accuracy. This 
                  option should be set to false if you receive warnings about 
                  negative charge densities after interpolation. 

    <b>Key name:</b>     freeze_occupied
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not to freeze the density and occupied 
                  orbitals after a restart. 

    <b>Key name:</b>     gw_residual_convergence_criterion
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Min value:</b>    1.000000e-14
    <b>Max value:</b>    4.000000e-04
    <b>Default:</b>      1.000000e-06
    <b>Description:</b>  The max value of the residual for unoccupied orbitals when 
                  performing a GW calculation. 

    <b>Key name:</b>     gw_residual_fraction
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.900000
    <b>Description:</b>  The residual value specified by gw_residual_convergence_criterion 
                  is applied to this fraction of the total spectrum. 

    <b>Key name:</b>     kohn_sham_ke_fft
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Special purpose flag which will force use of an FFT for the 
                  kinetic energy operator. 

    <b>Key name:</b>     kpoints
    <b>Required:</b>     yes
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  Normally kpoints are specified using the kpoint_mesh and 
                  kpoint_is_shift options but one can also enter a list of kpoints 
                  and their weights with this option. If kpoint_mesh is not 
                  specified or this is a bandstructure calculation this is required 
                  otherwise it is optional. 

    <b>Key name:</b>     kpoints_bandstructure
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  List of kpoints to use in a bandstructure calculation. For more 
                  detailed information look at the github wiki page on kpoint 
                  calculations. 

    <b>Key name:</b>     lcao_use_empty_orbitals
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Some pseudopotentials contain unbound atomic orbitals and this 
                  flag indicates whether or not they should be used for LCAO starts. 

    <b>Key name:</b>     md_steps_offset
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    0
    <b>Default:</b>      0
    <b>Description:</b>  

    <b>Key name:</b>     num_wanniers
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      0
    <b>Description:</b>  number of wannier functions to be used in wannier90 

    <b>Key name:</b>     rmg2bgw
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  Write wavefunction in G-space to BerkeleyGW WFN file. 

    <b>Key name:</b>     scf_steps_offset
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    0
    <b>Default:</b>      0
    <b>Description:</b>  

    <b>Key name:</b>     sqrt_interpolation
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not to use square root technique for 
                  density interpolation. 

    <b>Key name:</b>     total_scf_steps_offset
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    0
    <b>Default:</b>      0
    <b>Description:</b>  

    <b>Key name:</b>     use_cmix
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "true"
    <b>Description:</b>  Use adaptive interpolation can improve energies/forces for some 
                  but not all systems. 

    <b>Key name:</b>     use_cpdgemr2d
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  if set to true, we use Cpdgemr2d to change matrix distribution 

    <b>Key name:</b>     use_gpu_fd
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> Yes
    <b>Default:</b>      "false"
    <b>Description:</b>  Use gpus for kohn-sham orbital finite differencing. Depending on 
                  the balance of hardware characteristics this can provide a 
                  significant speedup but individual testing is required. 
                  Experimental. 

    <b>Key name:</b>     vxc_diag_nmax
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    10000
    <b>Default:</b>      1
    <b>Description:</b>  Maximum band index for diagonal Vxc matrix elements. 

    <b>Key name:</b>     vxc_diag_nmin
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    10000
    <b>Default:</b>      1
    <b>Description:</b>  Minimum band index for diagonal Vxc matrix elements. 

    <b>Key name:</b>     wannier90_scdm
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    -2147483647
    <b>Max value:</b>    2
    <b>Default:</b>      0
    <b>Description:</b>  use scdm method to set the trial wannier functions 

    <b>Key name:</b>     wannier90_scdm_sigma
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      1.000000
    <b>Description:</b>  when wannier90 is used to build wannier functions, the energy 
                  window parameter 

    <b>Key name:</b>     z_average_output_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "None"
    <b>Allowed:</b>      "wave functions" "potential and charge density" "None" 
    <b>Description:</b>  z_average_output_mode. 

</pre>
