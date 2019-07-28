
# How to compile
RMG builds have been tested using the GNU and PGI compilers as well as MKL.
Most development work is done using GNU which is the most reliable.
Cmake is used for configuration and out of source builds are preferred.
After cloning or downloading the repository the generic build instructions
are to first change into the top level directory and make a build subdir.

```Bash
cd rmgdft 
mkdir build
cd build
cmake ..
make -jN target
```

In this case N is an integer used for a parallel build and target specifies
the particular rmg module you wish to build. In the best case scenario things
will just work and on standard linux distributions this is often the case. On
non-standard clusters with complicated software stacks additional manual
configuration is often required. In particular one may have to specify a
particular set of modules to use and set some environment variables. For
example on the xsede machine comet as of July, 16, 2019 the following works.

```Bash
export CC=/opt/gnu/gcc/bin/gcc
export CXX=/opt/gnu/gcc/bin/c++
export FC=/opt/gnu/gcc/bin/gfortran
module load boost
module load fftw
module load cmake 
```

One can then execute
```Bash
cmake ..
make -jN target
```

Additional configuration may be performed at the cmake step. For example to
build with GPU support one should use.

```Bash
cmake -DRMG_GPU_ENABLED=1
```

Available targets include.
```
    rmg-cpu        Base code cpu only
    rmg-gpu        Base code with gpu support
    rmg-on-cpu     ON code cpu only
    rmg-on-gpu     ON code with gpu support
    rmg-negf-cpu   Non equilibrium greens function code cpu only
    rmg-negf-gpu   Non equilibrium greens function code with gpu support
    rmg-tddft-cpu  TDDFT based on rmg-on cpu only
    rmg-tddft-gpu  TDDFT based on rmg-on with gpu support
```
RMG GPU support at the present time is limited to Nvidia hardware and requires
Pascal or later hardware.


# Top level directory structure

cmake/
  Some additional cmake modules for finding specific libraries.

PlatformChecks/
  Code for checking whether certain features are suppored by compilers/libs.

SubprojectIncludes/
  Module specific cmake stuff.

RMG/
  Standard DFT code module (rmg-cpu and rmg-gpu binaries).

ON/
  Localized orbital DFT code module.

NEGF/
  Non equilibrium greens function code module.

TDDFT/
  Time dependent DFT code module.

RmgLib/
  Base code used by all modules.

US_PP/
  Code for working with pseudopotentials both US and NC.

Input/
  Routines for reading and parsing input files.

InternalPseudo/
  RMG includes a set of pseudopotentials built into the executable which
  are included here as compressed header files.

Misc/
  Miscellaneous code.

FiniteDiff/
  Higher level driver routines for finite differencing.

Force/
  Force routines.

MG/
  Multigrid routines. All multigrid functionality has been moved into RmgLib and this
  directory only has C bindings for the C++ class in RmgLib. Will eventually be deprecated.

Gpufuncs/
  Cuda code.

XC/
XC_useLIBXC/
  Interfaces for exchange correlation.

RMG_GUI/
  Gui setup code.

The next three directories contain 3rd party  libraries.

zfp/
  Compression library for floating point data.

scalapack/
  Parallel linear algegra/eigensolvers.

spglib/
  Symmetry routines.

Examples/
  Various examples.

Testing/
  Testing code. In development.
