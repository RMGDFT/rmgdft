# Internal pseudo potentials
RMG is normally built with internal pseudo potential support enabled which provides a set of pseudo potentials compiled into the executable file. These are the [GBRV 1.5](https://www.physics.rutgers.edu/gbrv/) ultrasoft pseudo potentials and are suitable for a wide range of problems. No additional action is required in order to use these.

# External pseudo potentials
Some types of calculations may require pseudo potentials with specific features not provided by the GBRV potentials. These include spin-orbit calculations, semi-local form and hybrid/exact exchange which the current version of RMG only implements for norm conserving potentials. These may be specified in the input file and RMG currently supports both the [UPF v2.0.1](http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format) and the [XML format](http://www.quantum-simulation.org/schemas/species.xsd).

# Semi-local pseudo potentials
As of version 4.1.0 RMG supports a semi-local pseudo potential implementation. Semi-local is useful for certain specific types of calculations but is normally much slower than a fully non-local implementation. If an input pseudopotential contains the necessary information to generate the a semi-local representation it can be enabled in RMG using by setting the input flag use_bessel_projectors="true".