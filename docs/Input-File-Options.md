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


## Control options

<pre>
    <b>Key name:</b>     a_length
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  First lattice constant. 

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
    <b>Allowed:</b>      "Exx Only" "NEB Relax" "Band Structure Only" "Psi Plot" "Plot" 
                  "Constant Pressure And Energy" "TDDFT" "Dimer Relax" "Constant 
                  Temperature And Energy" "Constant Volume And Energy" "Relax 
                  Structure" "Quench Electrons" 
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

    <b>Key name:</b>     description
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      ""
    <b>Allowed:</b>      
    <b>Description:</b>  Description of the run. 

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
    <b>Default:</b>      "Distributed fft"
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
    <b>Default:</b>      500
    <b>Description:</b>  Maximum number of self consistent steps to perform. Inner loop for 
                  hybrid functionals. 

    <b>Key name:</b>     noncollinear
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  if set true, calculate noncollinear 

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

    <b>Key name:</b>     spinorbit
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  if set true, calculate with spinorbit coupling 

    <b>Key name:</b>     start_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "LCAO Start"
    <b>Allowed:</b>      "Modified LCAO Start" "Restart TDDFT" "Start TDDFT" "Gaussian 
                  Start" "FIREBALL Start" "LCAO Start" "Restart From File" "Random 
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

    <b>Key name:</b>     vdw_corr
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "None"
    <b>Allowed:</b>      "DFT-D3" "DFT-D2" "Grimme-D2" "None" 
    <b>Description:</b>  Type of vdw correction 

    <b>Key name:</b>     vdwdf_grid_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Coarse"
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
    <b>Min value:</b>    5
    <b>Max value:</b>    50
    <b>Default:</b>      5
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
    <b>Allowed:</b>      "Tetragonal Primitive" "Cubic Body Centered" "Orthorhombic 
                  Primitive" "Cubic Face Centered" "Hexagonal Primitive" "Cubic 
                  Primitive" "None" 
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

    <b>Key name:</b>     lattice_units
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Bohr"
    <b>Allowed:</b>      "Angstrom" "Alat" "Bohr" 
    <b>Description:</b>  Units for the lattice vectors 

    <b>Key name:</b>     lattice_vector
    <b>Required:</b>     no
    <b>Key type:</b>     double array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "None"
    <b>Description:</b>  Lattice vectors, a0, a1, a2 

    <b>Key name:</b>     potential_grid_refinement
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    4
    <b>Default:</b>      0
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

    <b>Key name:</b>     use_symmetry
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  For non-gamma point, always true, for gamma point, optional 

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

    <b>Key name:</b>     filter_dpot
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  Flag indicating whether or not to filter density dependent 
                  potentials. 

    <b>Key name:</b>     filter_factor
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    0.060000
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.250000
    <b>Description:</b>  Filtering factor. 

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
    <b>Max value:</b>    10
    <b>Default:</b>      8
    <b>Description:</b>  RMG uses finite differencing to represent the kinetic energy 
                  operator and the accuracy of the representation is controllable by 
                  the kohn_sham_fd_order parameter. The default is 8 and is fine for 
                  most purposes but higher accuracy is obtainable with 10th order at 
                  the cost of some additional computational expense. 

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
    <b>Default:</b>      0.666667
    <b>Description:</b>  timestep for multigrid correction. 

    <b>Key name:</b>     kohn_sham_mucycles
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    6
    <b>Default:</b>      2
    <b>Description:</b>  Number of mu (also known as W) cycles to use in the kohn-sham 
                  multigrid preconditioner. 

    <b>Key name:</b>     kohn_sham_post_smoothing
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Min value:</b>    1
    <b>Max value:</b>    5
    <b>Default:</b>      2
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

</pre>
## Exchange correlation options

<pre>
    <b>Key name:</b>     exchange_correlation_type
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "AUTO_XC"
    <b>Allowed:</b>      "hartree-fock" "vdw-df-c09" "sla+pw+pbe+vdw1" "VDW-DF" "vdw-df" 
                  "gaupbe" "B3LYP" "hse" "mgga tb09" "AUTO_XC" "m06l" "VDW-DF-CX" 
                  "tpss" "ev93" "optbk88" "sogga" "wc" "HSE" "HCTH" "hcth" "Q2D" 
                  "q2d" "PBESOL" "tb09" "b86bpbe" "PW86PBE" "PBE0" "MGGA TB09" 
                  "pw86pbe" "REVPBE" "pbe" "revpbe" "GGA PBE" "BLYP" "pbe0" "pbesol" 
                  "blyp" "PBE" "GGA XP CP" "pw91" "GGA XB CP" "TB09" "optb86b" 
                  "olyp" "BP" "GGA BLYP" "bp" "b3lyp" "LDA" "vdw-df-cx" "PW91" "PZ" 
                  "pz" 
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
    <b>Description:</b>  order of Methefessel Paxton occupation. 

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
    <b>Allowed:</b>      "Error Function" "Gaussian" "Fermi Dirac" "MethfesselPaxton" "Cold 
                  Smearing" "Fixed" 
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

    <b>Key name:</b>     unoccupied_states_per_kpoint
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0
    <b>Max value:</b>    2147483647
    <b>Default:</b>      10
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
    <b>Max value:</b>    10
    <b>Default:</b>      5
    <b>Description:</b>  Number of previous steps to use when Broyden mixing is used to 
                  update the charge density. 

    <b>Key name:</b>     charge_broyden_scale
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.500000
    <b>Description:</b>  

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
    <b>Default:</b>      "Pulay"
    <b>Allowed:</b>      "Broyden" "Pulay" "Linear" 
    <b>Description:</b>  RMG supports Broyden, Pulay and Linear mixing When the davidson 
                  Kohn-Sham solver is selected Broyden or Pulay are preferred. For 
                  the multigrid solver Linear with potential acceleration is often 
                  (but not always) the best choice. 

    <b>Key name:</b>     charge_pulay_Gspace
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
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
    <b>Description:</b>  

    <b>Key name:</b>     charge_pulay_scale
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    1.000000
    <b>Default:</b>      0.500000
    <b>Description:</b>  

    <b>Key name:</b>     potential_acceleration_constant_step
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    4.000000
    <b>Default:</b>      0.000000e+00
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
    <b>Default:</b>      "Fast Relax"
    <b>Allowed:</b>      "LBFGS" "MD Min" "Quick Min" "FIRE" "Fast Relax" 
    <b>Description:</b>  Type of relaxation method to use for structural optimizations. 

    <b>Key name:</b>     renormalize_forces
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  Flag indicating whether or not to renormalize forces. 

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
    <b>Max value:</b>    512
    <b>Default:</b>      32
    <b>Description:</b>  Block size to use with scalapack. Optimal value is dependent on 
                  matrix size and system hardware. 

    <b>Key name:</b>     subdiag_driver
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "auto"
    <b>Allowed:</b>      "elpa" "cusolver" "auto" "scalapack" "magma" "lapack" 
    <b>Description:</b>  Driver type used for subspace diagonalization of the eigenvectors. 

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
    <b>Description:</b>  Block size to use when applying the non-local and S operators. 

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
    <b>Description:</b>  state_block used in nlforce. 

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
    <b>Description:</b>  Use asynchronous allreduce if available. 

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
    <b>Description:</b>  Use internal numa setup if available. 

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
                  Hubbard_U="Ni 6.5" 

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

    <b>Key name:</b>     tddft_mode
    <b>Required:</b>     no
    <b>Key type:</b>     string
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "electric field"
    <b>Allowed:</b>      "point charge" "electric field" 
    <b>Description:</b>  TDDFT mode 

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
    <b>Default:</b>      "Not done yet"
    <b>Description:</b>  cartesian coordinate of the point charge for tddft 

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
    <b>Default:</b>      0
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

</pre>
## Miscellaneous options

<pre>
    <b>Key name:</b>     E_POINTS
    <b>Required:</b>     no
    <b>Key type:</b>     integer
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    201
    <b>Max value:</b>    201
    <b>Default:</b>      201
    <b>Description:</b>  

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
    <b>Default:</b>      "true"
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

    <b>Key name:</b>     electric_field_magnitude
    <b>Required:</b>     no
    <b>Key type:</b>     double
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Min value:</b>    0.000000e+00
    <b>Max value:</b>    unlimited
    <b>Default:</b>      0.000000e+00
    <b>Description:</b>  Magnitude of external electric field. 

    <b>Key name:</b>     electric_field_vector
    <b>Required:</b>     no
    <b>Key type:</b>     double array
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "Not done yet"
    <b>Description:</b>  Components of the electric field. 

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
    <b>Description:</b>  if set true, calculate the exact exchange integrals 

    <b>Key name:</b>     fast_density
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
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

    <b>Key name:</b>     laplacian_autocoeff
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  if set to true, we use LaplacianCoeff.cpp to generate coeff 

    <b>Key name:</b>     laplacian_offdiag
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       Yes
    <b>Experimental:</b> No
    <b>Default:</b>      "false"
    <b>Description:</b>  if set to true, we use LaplacianCoeff.cpp to generate coeff 

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

    <b>Key name:</b>     use_cpdgemr2d
    <b>Required:</b>     no
    <b>Key type:</b>     boolean
    <b>Expert:</b>       No
    <b>Experimental:</b> No
    <b>Default:</b>      "true"
    <b>Description:</b>  if set to true, we use Cpdgemr2d to change matrix distribution 

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
    <b>Allowed:</b>      "potential and charge density" "wave functions" "None" 
    <b>Description:</b>  z_average_output_mode. 

</pre>
