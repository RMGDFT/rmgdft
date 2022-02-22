For the purposes of this document performance is defined as time to solution. Since RMG uses iterative methods time to solution depends on both the efficiency of the iterations and the speed with which they execute. These in turn may depend on both the hardware platform (e.g. cluster or workstation) and the problem type and size. There are many input file options and environment variables that may affect convergence rates and execution speed including.

## BLAS libaries
RMG is critically dependent on the performance of several double precision level 3 BLAS routines. In particular dgemm and dsyrk plus their complex equivalents zgemm and zsyrk for non-gamma calculations. Level 1 and 2 BLAS routines have an insignificant impact on RMG performance since they are rarely used. The matrix sizes passed to the routines will depend on the size of the problem (number of electronic wavefunctions and the real space basis size) as well as the number of MPI tasks used since RMG uses domain decomposition for the real space basis. As an example the initial electronic quench for a 512 atom NiO supercell run at an equivalent cutoff of 154 Rydbergs with 196 total MPI tasks is shown below.

| BLAS function | m    | n    | k     | Times called |
| ------------- | ----:| ----:| -----:| ------------:|
| DGEMM         | 3060  | 3570 | 31104| 329          |
| DGEMM         | 3060  | 3584 | 31104|   4          |
| DGEMM         | 3060  |   64 | 31104| 288          |
| DGEMM         | 31104 | 3570 | 3060 | 306          |
| DGEMM         | 31104 | 3570 | 3570 |  76          |
| DGEMM         | 31104 | 3584 | 3060 |   2          |
| DGEMM         | 31104 | 3584 | 3584 |   2          |
| DGEMM         | 3570  | 3570 | 31104|   4          |
| DSYRK         | --    | 3570 | 31104| 151          |
| DSYRK         | --    | 3584 | 31104|   4          |


## Charge density mixing
Pseudocode for the SCF cycle includes a mixing function <b>MIX</b> as shown below<br>

Initial charge density = <b>&rho;</b><sub>in</sub> and initial orbitals <b>&psi;</b><sub>i,j</sub> with i=0 and 1&lt;j&lt;N<br>
<b>do</b><br>
&nbsp;&nbsp;&nbsp;&nbsp;compute <b>V</b><sub>eff</sub>(<b>&rho;</b><sub>in</sub>)<br>
&nbsp;&nbsp;&nbsp;&nbsp;<b>for</b>(j=1 to N)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;multigrid preconditioner applied to <b>&psi;</b><sub>i,j</sub> using <b>V</b><sub>eff</sub>(<b>&rho;</b><sub>in</sub>)<br>
&nbsp;&nbsp;&nbsp;&nbsp;<b>end for</b><br>
&nbsp;&nbsp;&nbsp;&nbsp;diagonalize or orthogonalize <b>&psi;</b><sub>i,j</sub> => <b>&psi;</b><sub>i+1,j</sub><br>
&nbsp;&nbsp;&nbsp;&nbsp;generate <b>&rho;</b><sub>out</sub> from <b>&psi;</b><sub>i+1,j</sub><br>
&nbsp;&nbsp;&nbsp;&nbsp;update <b>&rho;</b><sub>in</sub> = <b>MIX</b>(<b>&rho;</b><sub>in</sub>, <b>&rho;</b><b><sub>out</sub>)<br>
&nbsp;&nbsp;&nbsp;&nbsp;<b>&Delta;&rho;</b><sub>i</sub> = <b>ABS</b>(<b>&rho;</b><sub>in</sub> - <b>&rho;</b><sub>out</sub>)<br>
&nbsp;&nbsp;&nbsp;&nbsp;i=i+1<br>
<b>while</b> (<b>&Delta;&rho;</b><sub>i</sub> &gt; tolerance)<br>

Three types of charge density mixing functions <b>MIX</b>(<b>&rho;</b><sub>in</sub>, <b>&rho;</b><sub>out</sub>) are available in RMG and are selected using the **charge_mixing_type** input option.

* **Linear** - This is the default and requires specifying a mixing constant <b>&alpha;</b>. New charge density is then calculated as <b>&rho;</b><sub>in</sub> = <b>&alpha;</b> <b>&rho;</b><sub>new</sub> + (1.0 - <b>&alpha;</b>) <b>&rho;</b><sub>old</sub>. Proper choice of <b>&alpha;</b> is crucial, for small values, reaching convergence may take many steps, while large values may lead to instability due to overshooting.
* **Pulay** - More sophisticated scheme, which uses charge densities from several previous steps to determine optimal charge density of the current step. may fail in some hard-to-converge cases. More info Pulay, Chem. Phys. Let. 73, 393 (1980).
* **Broyden** - Another multi-step scheme which is the preferred option when the Davidson Kohn-Sham solver is used.

## Potential acceleration
The SCF cycle outlined above recomputes <b>V</b><sub>eff</sub>(<b>&rho;</b><sub>in</sub>) once per SCF step (outer loop). It's possible to modify the cycle by updating <b>V</b><sub>eff</sub>(<b>&rho;</b><sub>in</sub>) in the inner loop over j. Specifically for each j we have a <b>&Delta;&psi;</b><sub>j</sub> = (<b>&psi;</b><sub>i+1,j</sub> - <b>&psi;</b><sub>i,j</sub>). Since the hartree potential <b>V</b><sub>h</sub> has a linear dependence on <<b>&psi;</b><sub>j</sub>|<b>&psi;</b><sub>j</sub>> we can compute an approximate update to <b>V</b><sub>h</sub> after computing each <b>&psi;</b><sub>i+1,j</sub> and then applying it to the preconditioning steps for successive <b>&psi;</b><sub>i+1,j</sub>. This has the effect of stabilizing the calculation and allows the use of a larger value of <b>&alpha;</b> when using linear charge density mixing (it is not compatible with Pulay or Broyden mixing). The input option **potential_acceleration_constant_step** is used to control this method and it's use is illustrated in the C60 examples.

## Subspace diagonalization driver
Depending on hardware resources and build configuration RMG users can select between several different diagonalizers. The optimal choice depends on the system hardware available and the problem size. Available options include.

* **Lapack** - Standard matrix algebra package. Required to build RMG. Is not parallel across MPI process's but can use multiple threads within a process via the BLAS libraries.
* **Scalapack** - Parallel version of lapack that decomposes matrices over MPI processes. A good choice for very large problems where the number of wavefunctions **N > 3000** but it may not be available for every hardware/software platform.
* **Cusolver** - GPU accelerated diagonalization routines. When suitable hardware is available this is often the best choice for problems where **N < 3000** but may not available for every hardware/software platform.


For small problems diagonalization normally comprises only a small part of the total execution time and the choice of driver is not critical. The computational work for diagonalization scales as **N<sup>3</sup>** though and as N increases this becomes a larger and larger portion of the calculation. With this in mind driver choice depends greatly on the hardware platform.

* **Workstation** - Lapack is only a good choice for workstations with a parallel BLAS library. If a parallel library is not available then Scalapack is preferred. Finally if a GPU is installed on the workstation Cusolver will usually be a better choice once **N** is more than a few hundred orbitals. 
* **Cluster** - The situation with clusters is considerably more complicated. The computational power provided by a single node has to be balanced against the communications speed available between nodes. Both Lapack and Cusolver are not parallel across nodes so their speed is limited by the speed of an individual node. Scalapack is parallel across nodes but is highly dependent on communications speed. Currently we have observed that for values of N < 3000 and both a high end communications fabric (Cray XK/XE series hardware) and high end GPU support (Nvidia Fermi/Kepler class hardware) Cusolver is usually the best choice. If GPU hardware is not available Scalapack is preferred. For **N > 3000** Scalapack will usually be the best choice for systems with fast communications.

* **Folded spectrum** A folded spectrum implementation is available that works in conjunction with the lapack or cusolver drivers. When enabled it can dramatically decrease the time required for subspace diagonalization on large systems. A paper describing the method may be downloaded from http://arxiv.org/abs/1502.07806 .