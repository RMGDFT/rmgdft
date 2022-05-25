# Internal pseudo potentials
RMG is normally built with internal pseudo potential support enabled which provides two sets of pseudo potentials compiled into the executable file. These are the [GBRV 1.5](https://www.physics.rutgers.edu/gbrv/) ultrasoft pseudo potentials and the [ONCVPSP](https://github.com/oncvpsp/oncvpsp) norm conserving pseudopotentials. ONCVPSP are the default and are suitable for a wide range of problems. No additional action is required in order to use them. Alternatively one can select the GBRV ultrasoft pseudo potentials by setting internal_pseudo_type = "ultrasoft" in the input file.

#  the External pseudo potentials
Some types of calculations may require pseudo potentials with specific features not provided by the included potentials. These include spin-orbit calculations and pseudopotentials in semi-local form. These may be specified in the input file and RMG currently supports both the [UPF v2.0.1](http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format) and the [XML format](http://www.quantum-simulation.org/schemas/species.xsd).

# Semi-local pseudo potentials
As of version 4.1.0 RMG supports a semi-local pseudo potential implementation. Semi-local is useful for certain specific types of calculations but is normally much slower than a fully non-local implementation. If an input pseudopotential contains the necessary information to generate the a semi-local representation it can be enabled in RMG using by setting the input flag use_bessel_projectors="true".
