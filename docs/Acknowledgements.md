The RMG source includes [spglib](https://atztogo.github.io/spglib/)[1] for working with crystal symmetries and the [libxc](https://tddft.org/programs/libxc/)[2] library of exchange correlation functionals. The [zfp](https://computing.llnl.gov/projects/floating-point-compression)[3] library is included and used for data compression/decompression while versions of the [LAPACK](http://www.netlib.org/lapack/) and [ScaLAPACK](http://www.netlib.org/scalapack/) libraries are included in order to help simplify the build process on certain platforms. A parallel FFT implementation from [LAMMPS](http://lammps.sandia.gov)[4] is also included as well as certain exchange-correlation routines from [Quantum Espresso](https://www.quantum-espresso.org/)[5].



1. Susi Lehtola, Conrad Steigemann, Micael J.T. Oliveira, and Miguel A.L. Marques, Recent developments in Libxc - A comprehensive library of functionals for density functional theory, Software X 7, 1 (2018)

2. “Spglib: a software library for crystal symmetry search”, Atsushi Togo and Isao Tanaka, https://arxiv.org/abs/1808.01590 (written at version 1.10.4)

3. Lindstrom, Peter. (2014). Fixed-Rate Compressed Floating-Point Arrays. IEEE Transactions on Visualization and Computer Graphics. 20. 10.1109/TVCG.2014.2346458. 

4. S. Plimpton, Fast Parallel Algorithms for Short-Range Molecular Dynamics, J Comp Phys, 117, 1-19 (1995)

5. P Giannozzi, O Andreussi, T Brumme, O Bunau, M Buongiorno Nardelli, M Calandra, R Car, C Cavazzoni, D Ceresoli, M Cococcioni, N Colonna, I Carnimeo, A Dal Corso, S de Gironcoli, P Delugas, R A DiStasio Jr, A Ferretti, A Floris, G Fratesi, G Fugallo, R Gebauer, U Gerstmann, F Giustino, T Gorni, J Jia, M Kawamura, H-Y Ko, A Kokalj, E Küçükbenli, M Lazzeri, M Marsili, N Marzari, F Mauri, N L Nguyen, H-V Nguyen, A Otero-de-la-Roza, L Paulatto, S Poncé, D Rocca, R Sabatini, B Santra, M Schlipf, A P Seitsonen, A Smogunov, I Timrov, T Thonhauser, P Umari, N Vast, X Wu and S Baroni, J.Phys.:Condens.Matter 29, 465901 (2017)