# History
RMG development began in 1993 and the original code was written in C with small portions in Fortran 77. Over time as features were added, this led to some issues with ease of maintenance so in 2013 the decision was made to improve the modularity of the code and rewrite significant portions in C++. The C++11 std (later updated to C++14) was targeted since it has native cross platform support for threading which is necessary for a hybrid MPI/threads implementation. This process is not entirely complete but as of July 2019 large portions of the code base have been rewritten in C++. (Update: Completed as of August 2019) Given constraints on manpower and time an incremental approach was chosen rather than a rewrite from scratch. Because of this, and inter operability requirements with modules still written mostly in C, some portions of the C++ code still have somewhat of a "C" flavor but this will diminish with future work.

Another development focus has been GPU integration. Currently RMG makes effective use of Nvidia GPU's to accelerate a number of computational tasks. Future work will provide support for AMD and Intel accelerators as well as extending the areas where GPU's are used in the code base.

# Toolchains
The primary development platform for RMG is GNU based. At various times it has been built with the Intel tool chain, the IBM toolchain and the PGI toolchain but GNU is the most reliable. RMG has supported native Windows builds in the past but with Microsoft's inclusion of Windows Subsystem for Linux in Windows 10 this has been deprecated. 

# Libraries
RMG performance is critically dependent on the availability of high performance BLAS libraries. On x86_64 platforms without GPU's both Openblas and Intel MKL work well. The generic libraries that ship with many linux distributions are quite slow and should not be used. When RMG is compiled with GPU support it will attempt to use the vendor supplied GPU libraries for BLAS operations.

RMG includes a number of libraries in the distribution. These include LAPACK, ScaLAPACK, spglib and the zfp libraries. Some of these may already be available on specific platforms but they are open source and are included to simplify the build process.

# RMG components
The build process can generate several different executable programs based on intended purpose and whether or not GPU support is available. These include.
* RMG Base which is a standard DFT code.
* RMG ON which is a localized orbital DFT variant.
* RMG NEGF which is Non Equilibrium Greens Function code.
* RMG TDDFT for Time Dependent DFT.

RMG Base is a production level code. The others are research/developer oriented codes that are less well documented and require significantly more code specific knowledge to use effectively.
