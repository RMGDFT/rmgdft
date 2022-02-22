# Environment variables
RMG is able to run in both pure MPI and hybrid mode. The mode is selected by setting environment variables that control the number of explicit threads used by RMG as well as implicit threading via OMP (both in RMG and in libraries that it links too such as BLAS).

Assuming that your shell is bash the number of RMG threads used is set with

`export RMG_NUM_THREADS=N ` 

The number of OMP threads is set with

`export OMP_NUM_THREADS=N`

These are typically set to the same value of N but in some cases better performance can be achieved by setting them to different values. When `N=1` RMG is running in pure MPI mode while a value greater than 1 selects hybrid mode. Setting either variable to a value greater than the number of cores available for each MPI process will lead to poor performance and should be avoided. The `OMP_WAIT_POLICY` environment variable may also impact performance and the best results are usually obtained using.

`export OMP_WAIT_POLICY=passive`
 
# Workstations
Modern workstations typically have 1-2 processors installed with each processor containing several physical CPU cores. With [hyperthreading](http://www.intel.com/content/www/us/en/architecture-and-technology/hyper-threading/hyper-threading-technology.html) or similar SMT technologies each physical core may appear as 2 or more logical cores to the OS and applications. For best performance the  number of MPI processes used during a calculation should not exceed the total number of logical cores. If the number of MPI processes used is less than the number of logical cores RMG will attempt to utilize the additional cores by using shared memory threads within each MPI process. The exact commands used to start a calculation vary depending on the specific MPI libraries your version of RMG was built with but a typical invocation for a version of RMG built using Openmpi would be.

mpirun -np 6 ./rmg-cpu in.c60

which specifies that 6 MPI processes will be started with input options contained in the file in.c60.

# Clusters
Unlike workstations, clusters are not commodities and their setup varies so widely that it's difficult to provide general running instructions. You should consult with your cluster administrators and review the documentation specific to your target platform for details on how to submit and run MPI jobs.

# Start mode
RMG offers several starting modes controlled by the start_mode input option. A full list of the available options may be found at the [Input File Options Page](https://github.com/RMGDFT/rmgdft/wiki/Input-File-Options) and can also be viewed from the command line using the --help option.

* **Random Start** - initial wavefunctions are set to random values. While useful in some situations an LCAO Start is usually a better choice.
* **LCAO Start** - initial wavefunctions are generated using a linear combination of atomic orbitals. Preferred starting mode.
* **Restart From File** - wavefunctions, charge density and all most other value are read from restart files.



# Calculation modes
Additional information about calculation modes can be found at the [Input File Options Page](https://github.com/RMGDFT/rmgdft/wiki/Input-File-Options).
* **Quench Electrons** - electron density and potentials are optimized variationally.
* **Relax Structures** - atomic structure is relaxed.
* **Band Structure Only** - band structures are calculated with fixed charge density and potentals. set kpoint_mesh="-1 1 1" and use the k point in the list of token kpoints="kx ky kz weight".
* **Constant Volume And Energy** - Molecular dynamics simulation with volume and energy held constant.

# Human readable output files

If the PLPLOT library is installed then image files in png format are created as well as the data files. The files are prefixed by the input file name so if that is in.c60 then.

**in.c60.rmsdv.xmgr**: convergence vs scf steps data
**in.c60.rmsdv.png**: convergence vs scf steps graph

**in.c60.bandstructure.xmgr**: Band structure data
**in.c60.bandstructure.png**: Band structure plot

**in.c60.dos.xmgr**: density of states data
**in.c60.dos.xmgr**: density of states plot
**in.c60.rho.xsf**:  file used to view atomic structure and charge density with XcrysDen.
**in.c60.log.\[nn\]**:  output file including energy, eigenvalue, force and convergence information. The value of nn starts from 00 and is incremented on successive runs up to 99.

The **rmsdv.png** file is updated dynamically during the run and can be monitored with an image viewer that refreshes automatically to watch the progress of the run. On a machine with ImageMagick installed the appropriate command would be.

**display -update 1 in.c60.rmsdv.png**

# Binary output (restart files)
RMG can generate binary output files suitable for restarts and interfacing with other codes in a variety of different formats. By default restarts and checkpoints are written in a parallel format where each MPI process writes it's own restart file. This is the fastest option but restart runs have to use the same number of MPI procs as the original run. A serial restart format is also available in which case RMG restarts can run using a different number of MPI procs than the original run as long as the number satisfies all other requirements for the run (e.g. for any given problem size there is an upper limit on the number of tasks that can be used).