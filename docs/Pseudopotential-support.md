# Internal pseudo potentials
RMG is normally built with internal pseudo potential support enabled which provides three sets of pseudo potentials compiled into the executable file. These are the [GBRV 1.5](https://www.physics.rutgers.edu/gbrv/) ultrasoft pseudo potentials, the [PseudoDojo](http://www.pseudo-dojo.org/) ONCV stringent pseudo potentials and the [SG15](http://www.quantum-simulation.org/potentials/sg15_oncv/) pseudo potentials generated using the [ONCVPP](http://www.mat-simresearch.com/) code. The SG15 variants are the default and are suitable for a wide range of problems. No additional action is required in order to use the SG15 potentials. Most of the SG15 potentials do not use core corrections and are compatible with hybrid XC functionals. The PseudDojo potentials do use core corrections and are incompatible with hybrid functionals. They may selected by adding the following line to your input file.<BR><BR>
`internal_pseudo_type="sg15"`<BR><BR>
`internal_pseudo_type="nc_accuracy"`<BR><BR>

 The GBRV potentials are softer than the SG15 or PseudoDojo potentials but are not suitable for use with hybrid functionals. They may be selected by adding the following line to your input file.<BR><BR>
`internal_pseudo_type="ultrasoft"`<BR><BR>

RMG also has support for using internally generated all-electron pseudopotentials using the method described here [All-Electron Plane-Wave Electronic Structure Calculations](https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c01191) by F. Gygi. This option is experimental and is intended to be used by experts/developers for pseudopotential verification. Unless you fall into that category you should not select this option which can be specified with.<BR><BR>
`internal_pseudo_type="all_electron"`<BR><BR>  

# External pseudo potentials
Some types of calculations may require pseudo potentials with specific features not provided by any of the internal potentials. These include spin-orbit calculations and some [pseudopotentials](https://pseudopotentiallibrary.org/) with many-body corrections. These may be specified in the input file and RMG currently supports both the [UPF v2.0.1](http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format) and the [XML format](http://www.quantum-simulation.org/schemas/species.xsd).

# Semi-local pseudo potentials
As of version 4.1.0 RMG supports a semi-local pseudo potential implementation. Semi-local is useful for certain specific types of calculations but is normally slower than a fully non-local implementation. If an input pseudopotential contains the necessary information to generate the a semi-local representation it can be enabled in RMG using by setting the input flag use_bessel_projectors="true".