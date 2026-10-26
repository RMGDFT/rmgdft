# QE and VASP input converters for RMG

Two command-line tools, also usable as Python modules, that write an RMG (rmg-dft) input file
from a Quantum ESPRESSO or a VASP input:

| file | what it does |
|---|---|
| `qe2rmg.py` | pw.x input (namelists and cards) to an RMG deck |
| `vasp2rmg.py` | POSCAR/CONTCAR, INCAR and KPOINTS to an RMG deck |
| `rmg_input.py` | shared code: deck model and writer, units, grid mapping, RMG lattice orientation, UPF header reader |
| `tests/` | pytest suite and small fixture inputs |

They need Python 3.8 or newer and numpy. The tests also need pytest. Neither ASE nor pymatgen is
used.

## Usage

```
qe2rmg.py pw.in -o input [--pseudo-dir DIR] [--grid NX NY NZ | --grid-spacing H] [--lenient]

vasp2rmg.py -d vasp_dir -o input [--pseudo Fe=Fe.upf --pseudo O=O.upf] [--potcar POTCAR]
            [--species "Fe O"] [--system-charge Q] [--grid NX NY NZ | --grid-spacing H] [--lenient]
```

- The deck goes to `-o` (stdout if omitted).
- RMG input files cannot contain comments. Its parser drops lines that start with `#` and
  truncates any line at `#`. A value also cannot contain `=` or `"`. So everything the converter
  changed, dropped or could not check goes to stderr and to `<output>.notes`. Read that file.
- The converter stops with exit status 2 when the input asks for physics that RMG, or this
  converter, cannot represent faithfully (DFT+U, electric fields, PAW UPF files, non-HSE06
  screening lengths, ...). `--lenient` drops such features instead, and lists each one in the
  notes as `DROPPED`.
- `--grid` sets the RMG `wavefunction_grid` directly. `--grid-spacing H` (in bohr) replaces the
  cutoff-based grid choice. `--grid-multiple N` (default 4) sets the divisor every grid
  dimension must have.
- `vasp2rmg.py` reads `POSCAR` (else `CONTCAR`), `INCAR`, `KPOINTS` and `POTCAR` from `-d` by
  default. Only the POSCAR is required.
- Spin-polarized RMG runs need an even number of MPI ranks (one group per spin channel).

Library use:

```python
from qe2rmg import convert
deck, notes = convert(open("pw.in").read(), input_dir=".")
open("input", "w").write(deck.render())
print(notes.text())
```

## How the main quantities are mapped

### Plane-wave cutoff to real-space grid

A grid of N points along a lattice vector of length L has spacing h = L/N. It resolves plane
waves up to |G| = pi/h, which is a kinetic-energy cutoff of (pi/h)^2 Ry. This is the
"Equivalent energy cutoff (psi)" that RMG prints in its log (`Misc/WriteHeader.cpp`). For each
lattice vector the converter therefore takes

    N_i >= L_i * sqrt(ecutwfc[Ry]) / pi        (i.e. h_i <= pi / sqrt(ecutwfc))

and rounds N_i up to the next multiple of 4 whose only prime factors are 2, 3, 5 and 7. RMG's
own automatic grid also uses multiples of 4, which keeps the multigrid levels usable; the
smoothness condition keeps the FFTs fast. RMG stops when the planar anisotropy of the grid is
above 1.5 (`Init.cpp`). If the rounded grid exceeds that limit, the converter refines the
coarsest direction until it does not, and says so in the notes.

VASP's `ENCUT` (eV) is mapped the same way, after dividing by 13.6057 eV/Ry. ENCUT is a PAW
cutoff, and RMG's norm-conserving pseudopotentials usually need a finer grid, so the notes
always warn about convergence. Without ENCUT, the largest POTCAR ENMAX is used (VASP's own
default), and failing that the converter asks for `--grid` or `--grid-spacing`.

The ratio QE `ecutrho / ecutwfc` becomes `potential_grid_refinement`. A fine grid r times
denser carries a density cutoff of r^2 times ecutwfc, so r = ceil(sqrt(ecutrho/ecutwfc)), with
a minimum of 2 (RMG's default, the equivalent of the norm-conserving dual of 4) and a maximum
of 4 (RMG's limit).

These choices were checked against the RMG binary. A converted Si deck and a converted HSE06
ZGNR deck were run with `rmg-cpu`. RMG's printed equivalent cutoffs (26.996 Ry and 63.677 Ry)
are the ones the converter predicts.

### Cell and atoms

The lattice is always written in bohr (`lattice_units = "Bohr"`, explicit `lattice_vector`). The
atoms are written in fractional coordinates (`atomic_coordinate_type = "Cell Relative"`).
Angstrom conversions use RMG's own constant (`A_a0 = 1.88972612499359`), so a deck in bohr means
exactly what RMG would build from the same input in angstrom.

RMG re-orients recognised Bravais lattices to its own convention (`Lattice::rotate_vectors`,
called from `Input/ReadCommon.cpp`). This keeps the lattice-vector lengths. It can, however,
change the angles between the vectors: QE's ibrav=3 body-centred cubic vectors become RMG's set,
with a different metric. It can also change the orientation of an orthorhombic or hexagonal cell
that was given in another orientation. As a result, neither the input fractional coordinates nor
the input Cartesian coordinates are always right in RMG's frame.

`rmg_input.orient_for_rmg` handles this. It holds a line-by-line transcription of RMG's
`lat2ibrav` and `rotate_vectors`. It finds the integer basis change M and the rotation Q with
`rotate(A) = M A Q`, and writes RMG's lattice with fractional coordinates `f M^-1`: the requested
crystal, rigidly rotated into RMG's frame. RMG's rotation is not idempotent for body-centred
cubic cells: its ibrav=3 output is re-detected as -3. The written cell is therefore iterated to
a fixed point, which RMG keeps unchanged. Explicit k-point lists are transformed with the same M.

Per-atom constraints (QE `if_pos`, VASP selective dynamics) become RMG's per-axis movable
columns. The initial spin goes into the RMG spin column. RMG takes a ratio r in [-0.5, 0.5], and
seeds the atom's density as (0.5 + r) for spin up and (0.5 - r) for spin down
(`US_PP/InitLocalObject.cpp`).

### Valence counts

QE `nbnd` and VASP `NBANDS` count all bands. RMG wants `unoccupied_states_per_kpoint` and adds
`ceil(nel/2)` occupied states (`Misc/init_states.cpp`, nel for noncollinear runs). That needs
the valence charge of every pseudopotential:

- `qe2rmg.py` reads `z_valence` from the UPF headers when the files are found.
- `vasp2rmg.py` reads it from the `--pseudo` UPF files.

When the charges are unknown, the band count is not converted, a warning is written, and RMG's
default applies.

## QE (pw.x) to RMG

| QE | RMG | notes |
|---|---|---|
| `calculation` scf / nscf / relax / vc-relax / md | `calculation_mode` "Quench Electrons" / "NSCF" / "Relax Structure" / "Relax Structure" + `cell_relax = "true"` / "Constant Volume And Energy" | md with `ion_temperature='nose'`: "Constant Temperature And Energy", Nose Hoover Chains, `md_temperature = tempw`. `bands` and `vc-md` are refused |
| `nstep` | `max_md_steps` | relax and md only (QE default 50) |
| `forc_conv_thr` (Ry/bohr) | `relax_max_force` (Ha/bohr) | divided by 2 |
| `ion_dynamics` bfgs / fire | `relax_method` LBFGS / FIRE | |
| `tstress` | `stress` | |
| `max_seconds` | `max_run_seconds` | key of the patch set (see below) |
| `ibrav` 1, 2, 3, -3, 4, 5, -5, 6, 7, 8, 9, -9, 10, 11, 12, -12, 14 + `celldm` or `A,B,C,cosAB,cosAC,cosBC` | `lattice_vector` (bohr) | QE `latgen` conventions; 13, -13 and 91 are refused (use ibrav=0) |
| `CELL_PARAMETERS` alat / bohr / angstrom | `lattice_vector` (bohr) | no option: alat if celldm(1) or A is set, else bohr |
| `ATOMIC_POSITIONS` alat / bohr / angstrom / crystal | `atoms`, Cell Relative | simple expressions such as `1/3` accepted; `crystal_sg` refused |
| `if_pos` | per-axis movable flags | |
| `ATOMIC_SPECIES` label, UPF | `pseudopotential` block, `pseudo_dir` | RMG reads UPF v2 directly; labels such as `Fe1` and `Fe2` map to element Fe; one UPF per element |
| `ecutwfc` | `wavefunction_grid` | see above; `--grid` and `--grid-spacing` override |
| `ecutrho` | `potential_grid_refinement` | |
| `input_dft` | `exchange_correlation_type` | pz/lda, pbe, pbesol, pw91, blyp, revpbe, tpss, m06l, tb09, pbe0, hse, b3lyp, gaupbe, hf, vdw-df family, scan/r2scan (libxc builds) |
| `exx_fraction` | `exx_fraction` | |
| `screening_parameter` | none | HSE: refused unless 0.106 bohr^-1, the value RMG always uses |
| `nqx1..3` | `exx_qpoint_mesh` | written for every hybrid (QE default 1 1 1; RMG's default is the full k mesh); must divide the k mesh |
| `x_gamma_extrapolation` | `x_gamma_extrapolation` | |
| `exxdiv_treatment` gygi-baldereschi / none | `exxdiv_treatment` | vcut_* refused |
| `exx_maxstep` | `max_exx_steps` | |
| `nspin = 2` | `spin_polarization = "true"` | |
| `starting_magnetization(i)` | per-atom spin column = m/2 | same up/down split of the atomic density |
| `noncolin`, `lspinorb`, `angle1`, `angle2` | `noncollinear`, `spinorbit`, 10-column atoms (r, theta, phi) | |
| `tot_magnetization` | `tot_magnetization` | smeared occupations; with fixed occupations: `states_count_and_occupation_spin_up/down` |
| `tot_charge` | `system_charge` | same sign convention (positive = electrons removed) |
| `occupations` smearing + `smearing` gaussian / m-p / m-v,cold / f-d | `occupations_type` MethfesselPaxton (`MP_order` 0) / MethfesselPaxton (`MP_order` 1) / "Cold Smearing" / "Fermi Dirac" | RMG's "Gaussian" occupation is a different function from QE's Gaussian, so MP order 0 (= 0.5 erfc) is used |
| `degauss` (Ry) | `occupation_electron_temperature_eV` | times 13.6057 |
| `occupations` tetrahedra / _lin / _opt | "Tetrahedron" + `tetra_method` Bloechl / Linear / Optimized | |
| `occupations = 'fixed'` | "Fixed" | |
| `nbnd` | `unoccupied_states_per_kpoint` | nbnd - ceil(nel/2); needs UPF z_valence |
| `K_POINTS automatic` n1 n2 n3 s1 s2 s3 | `kpoint_mesh`, `kpoint_is_shift` | same meaning: k = (i + s/2)/n |
| `K_POINTS gamma` | `kpoint_mesh "1 1 1"`, shift 0 | RMG then uses its gamma-point path |
| `K_POINTS tpiba / crystal` lists | `kpoints` (reciprocal-lattice units), `kpoint_mesh "-1 1 1"` | RMG normalises the weights and does not symmetrise the list |
| `conv_thr` (Ry) | `energy_convergence_criterion`, `rms_convergence_criterion`; hybrids also `exx_energy_convergence_criterion` | see "SCF and EXX stopping criteria" below |
| `electron_maxstep` | `max_scf_steps` | |
| `mixing_beta`, `mixing_ndim` | `charge_density_mixing`, `charge_broyden_order` | `charge_mixing_type = "Broyden"` (QE's plain mixing is Broyden) |
| `vdw_corr` grimme-d2 / grimme-d3, `dftd3_version`, `london` | `vdw_corr` DFT-D2 / DFT-D3, `dftd3_version` | |
| `nosym`, `noinv` | `use_symmetry = "0"`, `time_reversal = "false"` | |

### SCF and EXX stopping criteria

RMG stops an SCF loop as soon as either its energy test or its potential test is met
(`RMG/Common/Quench.cpp`, `Scf.cpp`). The energy test is the Hartree energy of the density residual,
1/2 ∫ (V_H[ρ_out] − V_H[ρ_in]) (ρ_out − ρ_in), plus for collinear spin the magnetization residual
(1/2π) ∫ |m_out − m_in|², the same measure QE's `conv_thr` bounds ("estimated scf accuracy",
`rho_ddot`). For noncollinear spin the energy test still sees only the total charge, so it can
report convergence while the magnetization is still relaxing. For a requested accuracy E (QE
`conv_thr`/2 or VASP `EDIFF`, in Ha) the converters therefore write

- without spin or with collinear spin: `energy_convergence_criterion` = E (between 1e-20 and RMG's
  maximum 1e-7), with `rms_convergence_criterion` = 1e-3 * sqrt(E in Ry) (between 1e-12 and 1e-4) as
  a backstop;
- with noncollinear spin: `energy_convergence_criterion` = 1e-20, RMG's floor, so the RMS criterion
  above decides;
- for hybrids, `exx_energy_convergence_criterion` = E: the EXX loop stops when
  |f1 - (f0+f2)/2| < E, the quantity QE calls dexx and compares with conv_thr.

RMG builds without the magnetization term (before this pull request) stop collinear spin-polarized
runs too early: a 2-ZGNR PBE run stopped at 4e-13 Ha with RMS[dV] = 2e-4 and a total moment of 1.00016
where 1 was imposed, and spin-polarized hybrid runs stalled at |exx_delta| ~ 1e-5 Ha. With such a
build, set `energy_convergence_criterion` to 1e-20 by hand for spin-polarized decks.

## VASP to RMG

| VASP | RMG | notes |
|---|---|---|
| POSCAR (VASP 5 species line, or VASP 4 + `--species`/POTCAR), Direct/Cartesian | `lattice_vector`, `atoms` (Cell Relative) | species such as `Fe_pv` map to Fe |
| scale factor, negative scale (= volume), three scale factors | applied to lattice and Cartesian positions | |
| Selective dynamics | per-axis movable flags | |
| POTCAR | none | PAW cannot be converted; read only for TITEL (species), ZVAL (NELECT) and ENMAX (default ENCUT) |
| `--pseudo El=file.upf` | `pseudopotential`, `pseudo_dir` | every element or none; without it RMG uses its internal norm-conserving pseudopotentials |
| `ENCUT` (eV) | `wavefunction_grid` | see above |
| `ISPIN = 2` | `spin_polarization` | |
| `MAGMOM` (incl. `N*value`) | per-atom spin column = MAGMOM / (2 z_valence) | without UPF valence: MAGMOM/10, clamped to 0.5, a starting guess only |
| `LNONCOLLINEAR`, `LSORBIT`, 3-component `MAGMOM` | `noncollinear`, `spinorbit`, (r, theta, phi) | `SAXIS` other than 0 0 1 refused |
| `NUPDOWN` | `tot_magnetization` | |
| `GGA` PE / PS / 91 / RE / CA | `exchange_correlation_type` pbe / pbesol / pw91 / revpbe / pz | no GGA: functional from the pseudopotentials (AUTO_XC) |
| `METAGGA` SCAN / R2SCAN / TPSS / M06L / MBJ | scan / r2scan / tpss / m06l / tb09 | |
| `LHFCALC` (+ `HFSCREEN = 0.2`, + `GGA = B3`) | pbe0 (hse, b3lyp) | HFSCREEN 0.2 A^-1 = 0.1058 bohr^-1, taken as RMG's HSE06 (0.106); other screening lengths and custom ALDAC/AGGAX/AGGAC refused |
| `AEXX` | `exx_fraction` | when it differs from the functional's default |
| `NKRED`, `NKREDX/Y/Z` | `exx_qpoint_mesh` = k mesh / NKRED | |
| `ISMEAR` -5 / -4 / -1 / 0 / N>0 | "Tetrahedron" Bloechl / Linear, "Fermi Dirac", MethfesselPaxton `MP_order` 0 / N | VASP defaults ISMEAR 1 and SIGMA 0.2 are written explicitly |
| `SIGMA` (eV) | `occupation_electron_temperature_eV` | |
| `NBANDS` | `unoccupied_states_per_kpoint` | needs `--pseudo` (valence) |
| `NELECT` | `system_charge` = sum(POTCAR ZVAL) - NELECT | needs POTCAR, or `--system-charge` |
| `EDIFF` (eV) | `energy_convergence_criterion`, `rms_convergence_criterion`; hybrids also `exx_energy_convergence_criterion` | see "SCF and EXX stopping criteria" below |
| `NELM` | `max_scf_steps` | |
| `AMIX` | `charge_density_mixing` | |
| `IBRION` 1/2/3 + `NSW` | "Relax Structure", `max_md_steps` | `ISIF = 3`: `cell_relax`; ISIF 4-7 refused |
| `EDIFFG < 0` (eV/A) | `relax_max_force` (Ha/bohr) | |
| `IBRION = 0` + `SMASS` -3 / >= 0 (`MDALGO` 0 or 2) | "Constant Volume And Energy" / "Constant Temperature And Energy" with `md_temperature = TEBEG` | |
| `ICHARG = 11` | "NSCF" | |
| `IVDW` 1,10 / 11 / 12 | `vdw_corr` DFT-D2 / DFT-D3 (`dftd3_version` 3) / DFT-D3 (`dftd3_version` 4) | |
| `ISYM` 0 / -1 | `use_symmetry = "0"` (+ `time_reversal = "false"`) | |
| KPOINTS Gamma / Monkhorst-Pack + shifts | `kpoint_mesh`, `kpoint_is_shift` | MP grids are shifted by half a step for even N; total shifts must be 0 or 1/2 |
| KPOINTS fully automatic (`Auto`, R_k) | `kpoint_mesh` | N_i = max(1, int(R_k b_i + 0.5)), Gamma-centred |
| no KPOINTS: `KSPACING`, `KGAMMA` | `kpoint_mesh`, `kpoint_is_shift` | N_i = max(1, ceil(2 pi b_i / KSPACING)); VASP default 0.5 |
| KPOINTS explicit list (reciprocal) | `kpoints` | |

## What is not converted, and why

These features stop the conversion (or are dropped with `--lenient`):

- **DFT+U** (QE `lda_plus_u`/`HUBBARD`, VASP `LDAU*`). RMG's `Hubbard_U` needs orbital labels and
  uses its own projectors, so a literal translation would not reproduce the same U.
- **Electric fields, Berry phase, dipole and gate corrections, ESM, FCP, RISM, `assume_isolated`,
  VASPsol, linear response** (`LEPSILON`, `LOPTICS`, `LCALCEPS`, ...). There is no RMG counterpart
  in this converter.
- **Screening lengths other than HSE06's.** RMG has no input for the HSE screening parameter; it
  always uses 0.106 bohr^-1.
- **`exxdiv_treatment` vcut_ws / vcut_spherical.** RMG offers only gygi-baldereschi and none.
- **Custom hybrid mixing** (VASP `ALDAC`/`AGGAX`/`AGGAC`, meta-GGA hybrids).
- **PAW and UPF v1 files.** RMG reads norm-conserving and ultrasoft UPF v2 (XML) only. Convert UPF
  v1 with QE's `upfconv.x -u`. VASP POTCARs are never converted.
- **Band-structure paths** (QE `tpiba_b`/`crystal_b`, VASP line mode) and Cartesian VASP k-lists.
  Write an RMG `kpoints_bandstructure` block by hand.
- **Hybrid functionals with explicit k-point lists.** RMG's exact exchange needs a k mesh.
- **k-point shifts other than 0 or 1/2 of a grid step.**
- **QE ibrav 13, -13, 91 and `crystal_sg` positions.** Use ibrav=0.
- **Body-centred tetragonal cells that RMG detects as such.** RMG's `rotate_vectors` re-orients
  them using the primitive-vector length for both a and c, which is a different lattice.
- **Constrained cell shapes** (QE `cell_dofree` other than all, VASP ISIF 4-7), **target
  pressure, temperature ramps, VASP `MDALGO` other than 0/2, `SMASS = -1`, phonon modes**
  (`IBRION` 5-8, 44), **constrained magnetization, `CONSTRAINTS`, `OCCUPATIONS` card
  (`from_input`), ISMEAR -2/-3**.

These are left out with a note or warning:

- **Run control** (prefix, outdir, verbosity, NPAR/KPAR/NCORE, LWAVE, PREC, ALGO, ...). RMG has its
  own I/O and parallelisation settings.
- **MD time step** (QE `dt`, VASP `POTIM`). Set `ionic_time_step` by hand.
- **Atomic masses** from ATOMIC_SPECIES. RMG uses its own table.
- **QE `ecutfock`, `ecutvcut`, `ace`, `dftd3_threebody`**, and VASP `PRECFOCK`, `GGA_COMPAT`.
- **The VASP treatment of the Coulomb singularity in exact exchange.** RMG's defaults
  (gygi-baldereschi, x_gamma_extrapolation) apply.
- **QE `etot_conv_thr` and VASP `EDIFFG > 0`.** RMG stops relaxations on forces only.
- **Unknown keys and tags.** A warning names each one.

## Recent RMG keys

`tot_magnetization`, `exx_qpoint_mesh`, `max_run_seconds` and `exx_energy_convergence_criterion`
are registered in `Input/ReadCommon.cpp` by the same change set as these converters. RMG ignores
unknown keys silently, so an older build would run with a free moment, the full k mesh for exact
exchange, no time limit and only the vexx_rms EXX stop. Whenever the converter writes one of these keys, it says so in the
notes.

## Tests

```
cd scripts/converters
python3 -m pytest -q
```

The suite covers:

- every unit and format branch: alat, bohr, angstrom and crystal positions; CELL_PARAMETERS
  units; celldm against A,B,C; ibrav lattices; POSCAR scaling (including negative and
  per-axis scales); selective dynamics;
- k-point meshes and shifts;
- MAGMOM expansion;
- the functional, hybrid and occupation mappings;
- the grid mapping, checked against an RMG production deck;
- RMG's lattice re-orientation;
- the UPF reader.

Round-trip tests read the written deck back with `rmg_input.parse_rmg_deck` and compare the cell
and the Cartesian positions with the source input to 1e-8 bohr. Refusal tests check that
unsupported input stops with `ConversionError`, or is dropped under `--lenient`.

The UPF files in `tests/data/pseudo` are header-only fixtures and cannot be used for
calculations.
