"""Tests for the Quantum ESPRESSO pw.x -> RMG converter."""

import math
import os

import numpy as np
import pytest

from conftest import DATA, PSEUDO, read_data
from qe2rmg import convert, eval_number, main, parse_pw_input, qe_lattice
from rmg_input import (BOHR_PER_ANGSTROM, RY_EV, ConversionError, deck_structure,
                       parse_rmg_deck)


def run(text, **kw):
    kw.setdefault("input_dir", DATA)
    deck, notes = convert(text, **kw)
    return parse_rmg_deck(deck.render()), notes


def mini(system="", control="", electrons="", positions=None, cell=None, kpoints=None,
         species=None, ibrav=1, nat=2):
    """A small Si input with replaceable parts (pseudo_dir points at the test UPF stubs)."""
    sysblock = "  ibrav = %d\n" % ibrav
    if ibrav == 1:
        sysblock += "  celldm(1) = 10.0\n"
    text = ("&CONTROL\n  pseudo_dir = 'pseudo'\n%s/\n&SYSTEM\n%s  nat = %d, ntyp = 1\n"
            "  ecutwfc = 30.0\n%s/\n&ELECTRONS\n%s/\n" % (control, sysblock, nat, system,
                                                         electrons))
    text += species or "ATOMIC_SPECIES\n  Si 28.0855 Si.header-only.UPF\n"
    if cell:
        text += cell
    text += positions or "ATOMIC_POSITIONS crystal\n  Si 0 0 0\n  Si 0.25 0.25 0.25\n"
    text += kpoints or "K_POINTS automatic\n  2 2 2 0 0 0\n"
    return text


# ---------------------------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------------------------


def test_namelist_parsing_details():
    text = """
&control
   calculation='relax', pseudo_dir = './a/b!c/'   ! comment with = and /
   tstress = .TRUE. , tprnfor=T
/
&SYSTEM
  ibrav=2, celldm(1)=1.02d1, nat=2,ntyp=1
  starting_magnetization(1) = 0.5, 0.25
  ecutwfc = 3.0D+01
/
ATOMIC_SPECIES
  Si 28.0 Si.UPF
"""
    pw = parse_pw_input(text)
    c, s = pw.nl("CONTROL"), pw.nl("SYSTEM")
    assert c["calculation"] == "relax"
    assert c["pseudo_dir"] == "./a/b!c/"
    assert c["tstress"] is True and c["tprnfor"] is True
    assert s["celldm(1)"] == 10.2 and s["nat"] == 2 and s["ecutwfc"] == 30.0
    assert s["starting_magnetization(1)"] == 0.5 and s["starting_magnetization(2)"] == 0.25
    assert pw.cards["ATOMIC_SPECIES"].lines == ["Si 28.0 Si.UPF"]


@pytest.mark.parametrize("txt,val", [("1/3", 1 / 3), ("-0.5d0", -0.5), ("2^2", 4.0),
                                     ("(1+1)/4", 0.5), ("1.0E-2", 0.01)])
def test_eval_number(txt, val):
    assert eval_number(txt) == pytest.approx(val)


def test_eval_number_rejects_code():
    with pytest.raises(ConversionError):
        eval_number("__import__('os')")


# ---------------------------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------------------------


def test_qe_lattices_common_ibrav():
    a = 10.0
    assert np.allclose(qe_lattice(1, {1: a}), np.eye(3) * a)
    assert np.allclose(qe_lattice(2, {1: a}), np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]]) * a / 2)
    assert np.allclose(qe_lattice(3, {1: a}), np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]]) * a / 2)
    assert np.allclose(qe_lattice(4, {1: a, 3: 1.6}),
                       [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, 1.6 * a]])
    assert np.allclose(qe_lattice(6, {1: a, 3: 1.3}), np.diag([a, a, 1.3 * a]))
    assert np.allclose(qe_lattice(8, {1: a, 2: 1.2, 3: 1.5}), np.diag([a, 1.2 * a, 1.5 * a]))


def test_qe_lattice_volumes_and_shapes():
    a = 7.0
    # fcc, bcc, trigonal (5 and -5 are the same lattice), triclinic
    assert abs(np.linalg.det(qe_lattice(2, {1: a}))) == pytest.approx(a ** 3 / 4)
    assert abs(np.linalg.det(qe_lattice(3, {1: a}))) == pytest.approx(a ** 3 / 2)
    t5, tm5 = qe_lattice(5, {1: a, 4: 0.3}), qe_lattice(-5, {1: a, 4: 0.3})
    assert np.allclose(t5 @ t5.T, tm5 @ tm5.T)
    assert np.allclose(np.linalg.norm(t5, axis=1), a)
    cal, cbe, cga = 0.1, 0.2, 0.3
    tri = qe_lattice(14, {1: a, 2: 1.1, 3: 1.2, 4: cal, 5: cbe, 6: cga})
    l = np.linalg.norm(tri, axis=1)
    assert l == pytest.approx([a, 1.1 * a, 1.2 * a])
    assert tri[1] @ tri[2] / l[1] / l[2] == pytest.approx(cal)
    assert tri[0] @ tri[2] / l[0] / l[2] == pytest.approx(cbe)
    assert tri[0] @ tri[1] / l[0] / l[1] == pytest.approx(cga)


def test_celldm_and_abc_give_the_same_deck():
    a_ang = 5.0
    cd = a_ang * BOHR_PER_ANGSTROM
    t1 = mini(system="  celldm(2) = 1.2, celldm(3) = 1.5\n", ibrav=8).replace(
        "  ibrav = 8\n", "  ibrav = 8\n  celldm(1) = %.15f\n" % cd)
    t2 = mini(system="  A = 5.0, B = 6.0, C = 7.5\n", ibrav=8)
    p1, _ = run(t1)
    p2, _ = run(t2)
    assert np.allclose(deck_structure(p1).lattice_bohr, deck_structure(p2).lattice_bohr, atol=1e-8)
    assert np.allclose(deck_structure(p2).lattice_bohr, np.diag([5.0, 6.0, 7.5]) * BOHR_PER_ANGSTROM,
                       atol=1e-8)


def test_celldm_and_abc_together_are_refused():
    with pytest.raises(ConversionError):
        run(mini(system="  A = 5.0\n"))


@pytest.mark.parametrize("ibrav", [13, -13, 91])
def test_unsupported_ibrav_is_refused(ibrav):
    with pytest.raises(ConversionError):
        run(mini(ibrav=ibrav, system="  celldm(1) = 10.0, celldm(2) = 1.1, celldm(3) = 1.2\n"))


def test_cell_parameters_units_agree():
    pos = "ATOMIC_POSITIONS crystal\n  Si 0 0 0\n  Si 0.3 0.2 0.1\n"
    a = 4.0
    ang = mini(ibrav=0, positions=pos,
               cell="CELL_PARAMETERS angstrom\n %f 0 0\n 0.5 %f 0\n 0.2 0.3 %f\n" % (a, a, a))
    bohr = mini(ibrav=0, positions=pos,
                cell="CELL_PARAMETERS bohr\n %.15f 0 0\n %.15f %.15f 0\n %.15f %.15f %.15f\n"
                % tuple(x * BOHR_PER_ANGSTROM for x in (a, 0.5, a, 0.2, 0.3, a)))
    alat = mini(ibrav=0, positions=pos, system="  celldm(1) = %.15f\n" % (a * BOHR_PER_ANGSTROM),
                cell="CELL_PARAMETERS alat\n 1 0 0\n %.15f 1 0\n %.15f %.15f 1\n"
                % (0.5 / a, 0.2 / a, 0.3 / a))
    structs = [deck_structure(run(t)[0]) for t in (ang, bohr, alat)]
    for s in structs[1:]:
        assert np.allclose(s.lattice_bohr, structs[0].lattice_bohr, atol=1e-8)
        assert np.allclose(s.cart_bohr, structs[0].cart_bohr, atol=1e-8)


def test_position_units_agree():
    """The same atoms given in crystal, alat, bohr and angstrom units give one deck."""
    a = 10.0
    frac = np.array([[0.0, 0.0, 0.0], [0.27, 0.31, 0.22]])
    cart = frac * a
    variants = {
        "crystal": frac,
        "alat": cart / a,
        "bohr": cart,
        "angstrom": cart / BOHR_PER_ANGSTROM,
    }
    ref = None
    for opt, xyz in variants.items():
        pos = "ATOMIC_POSITIONS {%s}\n" % opt + "".join(
            "  Si %.15f %.15f %.15f\n" % tuple(r) for r in xyz)
        st = deck_structure(run(mini(positions=pos))[0])
        assert np.allclose(st.frac, frac, atol=1e-10), opt
        if ref is None:
            ref = st.cart_bohr
        assert np.allclose(st.cart_bohr, ref, atol=1e-8)


def test_position_expressions_and_if_pos():
    pos = "ATOMIC_POSITIONS crystal\n  Si 0 0 0 0 0 0\n  Si 1/4 1/4 1/4 1 0 1\n"
    p, notes = run(mini(positions=pos, control="  calculation = 'relax'\n"))
    rows = [l.split() for l in p["atoms"].splitlines()]
    assert rows[0][4:7] == ["0", "0", "0"] and rows[1][4:7] == ["1", "0", "1"]
    assert float(rows[1][1]) == 0.25
    assert any("if_pos" in w for w in notes.warnings)


# ---------------------------------------------------------------------------------------------
# Fixtures: round trips
# ---------------------------------------------------------------------------------------------


def test_si_diamond_fixture():
    p, notes = run(read_data("si_diamond.pw.in"))
    st = deck_structure(p)
    a = 10.26
    qe = np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]]) * a / 2
    # RMG re-orients fcc; the metric (lengths and angles) and the fractions are unchanged
    assert np.allclose(st.lattice_bohr @ st.lattice_bohr.T, qe @ qe.T, atol=1e-8)
    assert np.allclose(st.frac, [[0, 0, 0], [0.25, 0.25, 0.25]], atol=1e-10)
    assert np.linalg.norm(st.cart_bohr[1] - st.cart_bohr[0]) == pytest.approx(a * math.sqrt(3) / 4)
    assert p["calculation_mode"] == "Quench Electrons"
    assert p["kpoint_mesh"] == "4 4 4" and p["kpoint_is_shift"] == "1 1 1"
    assert p["unoccupied_states_per_kpoint"] == "4"     # nbnd 8 - 8 electrons / 2
    assert p["occupations_type"] == "Fixed"
    assert p["wavefunction_grid"] == "16 16 16"
    assert p["potential_grid_refinement"] == "2"
    # no spin: conv_thr 1e-10 Ry maps onto RMG's own energy test (the same Hartree-residual measure), in Ha
    assert float(p["energy_convergence_criterion"]) == pytest.approx(5e-11)
    assert float(p["rms_convergence_criterion"]) == pytest.approx(1e-8)
    assert p["pseudopotential"].split() == ["Si", "Si.header-only.UPF"]
    assert "exchange_correlation_type" not in p and "spin_polarization" not in p
    assert not notes.warnings


def test_zgnr_hse_fixture():
    text = read_data("zgnr_hse.pw.in")
    p, notes = run(text)
    st = deck_structure(p)
    pw = parse_pw_input(text)
    cell = np.array([[float(x) for x in l.split()] for l in pw.cards["CELL_PARAMETERS"].lines])
    pos = np.array([[float(x) for x in l.split()[1:4]] for l in pw.cards["ATOMIC_POSITIONS"].lines])
    assert np.allclose(st.lattice_bohr, cell * BOHR_PER_ANGSTROM, atol=1e-8)
    assert np.allclose(st.cart_bohr, pos * BOHR_PER_ANGSTROM, atol=1e-8)
    assert st.symbols == ["C"] * 6 + ["H"] * 2
    assert [m[0] for m in st.moments] == [0.25] * 6 + [0.0] * 2
    assert p["spin_polarization"] == "true"
    assert p["exchange_correlation_type"] == "hse"
    assert p["exx_qpoint_mesh"] == "1 1 4"
    assert p["x_gamma_extrapolation"] == "true"
    assert p["exxdiv_treatment"] == "gygi-baldereschi"
    assert p["kpoint_mesh"] == "1 1 8" and p["kpoint_is_shift"] == "0 0 1"
    assert float(p["tot_magnetization"]) == 0.5
    assert p["occupations_type"] == "Cold Smearing"
    assert float(p["occupation_electron_temperature_eV"]) == pytest.approx(0.002 * RY_EV)
    assert p["unoccupied_states_per_kpoint"] == "7"     # nbnd 20 - (6*4 + 2*1) / 2
    assert p["potential_grid_refinement"] == "2"
    assert p["max_scf_steps"] == "300"
    # conv_thr 1e-9 Ry: EXX stops on |exx_delta| < 5e-10 Ha; spin-polarized, so the energy test
    # (blind to the magnetization residual) is at RMG's floor and RMS[dV] < 1e-3*sqrt(1e-9) stops the SCF
    assert float(p["exx_energy_convergence_criterion"]) == pytest.approx(5e-10)
    assert float(p["energy_convergence_criterion"]) == pytest.approx(1e-20)
    assert float(p["rms_convergence_criterion"]) == pytest.approx(1e-3 * math.sqrt(1e-9))
    assert float(p["charge_density_mixing"]) == 0.3
    # grid: every spacing at or below pi/sqrt(60 Ry)
    g = [int(x) for x in p["wavefunction_grid"].split()]
    h = np.linalg.norm(st.lattice_bohr, axis=1) / g
    assert np.all(h <= math.pi / math.sqrt(60.0) + 1e-12)
    assert g == [96, 72, 12]
    assert "3-ZGNR HSE06 FM test" in p["description"]
    assert not notes.warnings


# ---------------------------------------------------------------------------------------------
# Functionals, exact exchange, occupations
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("dft,xc,hyb", [("PBE", "pbe", False), ("pz", "pz", False),
                                        ("lda", "pz", False), ("pbesol", "pbesol", False),
                                        ("hse", "hse", True), ("PBE0", "pbe0", True),
                                        ("b3lyp", "b3lyp", True), ("hf", "hartree-fock", True)])
def test_functional_mapping(dft, xc, hyb):
    p, _ = run(mini(system="  input_dft = '%s'\n" % dft))
    assert p["exchange_correlation_type"] == xc
    assert ("exx_qpoint_mesh" in p) == hyb
    if hyb:
        assert p["exx_qpoint_mesh"] == "1 1 1"           # QE default nqx
        assert p["x_gamma_extrapolation"] == "true"


def test_unknown_functional_is_refused():
    with pytest.raises(ConversionError):
        run(mini(system="  input_dft = 'm06-2x-nonsense'\n"))


def test_hybrid_qmesh_must_divide_kmesh():
    with pytest.raises(ConversionError):
        run(mini(system="  input_dft = 'pbe0', nqx1 = 3\n"))
    p, _ = run(mini(system="  input_dft = 'pbe0', nqx1 = 2, nqx2 = 1, nqx3 = 2\n"))
    assert p["exx_qpoint_mesh"] == "2 1 2"


def test_hse_screening_and_divergence_options():
    with pytest.raises(ConversionError):
        run(mini(system="  input_dft = 'hse', screening_parameter = 0.2\n"))
    with pytest.raises(ConversionError):
        run(mini(system="  input_dft = 'hse', exxdiv_treatment = 'vcut_ws'\n"))
    p, notes = run(mini(system="  input_dft = 'hse', exxdiv_treatment = 'vcut_ws'\n"),
                   lenient=True)
    assert "exxdiv_treatment" not in p and any("DROPPED" in w for w in notes.warnings)
    p, _ = run(mini(system="  input_dft = 'hse', x_gamma_extrapolation = .false., "
                           "exx_fraction = 0.3\n  exxdiv_treatment = 'none'\n",
                    electrons="  exx_maxstep = 7\n"))
    assert p["x_gamma_extrapolation"] == "false" and p["exxdiv_treatment"] == "none"
    assert float(p["exx_fraction"]) == 0.3 and p["max_exx_steps"] == "7"


@pytest.mark.parametrize("smearing,otype,order", [
    ("gaussian", "MethfesselPaxton", "0"), ("m-p", "MethfesselPaxton", "1"),
    ("marzari-vanderbilt", "Cold Smearing", None), ("fermi-dirac", "Fermi Dirac", None)])
def test_smearing_mapping(smearing, otype, order):
    p, _ = run(mini(system="  occupations = 'smearing', smearing = '%s', degauss = 0.01\n"
                           % smearing))
    assert p["occupations_type"] == otype
    assert p.get("MP_order") == order
    assert float(p["occupation_electron_temperature_eV"]) == pytest.approx(0.01 * RY_EV)


@pytest.mark.parametrize("occ,method", [("tetrahedra", "Bloechl"), ("tetrahedra_lin", "Linear"),
                                        ("tetrahedra_opt", "Optimized")])
def test_tetrahedra(occ, method):
    p, _ = run(mini(system="  occupations = '%s'\n" % occ))
    assert p["occupations_type"] == "Tetrahedron" and p["tetra_method"] == method


def test_spin_fixed_occupations():
    with pytest.raises(ConversionError):
        run(mini(system="  nspin = 2, starting_magnetization(1) = 0.5\n"))
    p, _ = run(mini(system="  nspin = 2, tot_magnetization = 2, nbnd = 7\n"))
    # 8 electrons, M = 2: 5 up, 3 down, 7 states per spin
    assert p["states_count_and_occupation_spin_up"] == "5 1.0 2 0.0"
    assert p["states_count_and_occupation_spin_down"] == "3 1.0 4 0.0"
    assert "tot_magnetization" not in p


def test_noncollinear_spin_orbit():
    p, notes = run(mini(system="  noncolin = .true., lspinorb = .true.\n"
                               "  starting_magnetization(1) = 0.4, angle1(1) = 90, angle2(1) = 45\n"
                               "  occupations = 'smearing', degauss = 0.01\n"))
    assert p["noncollinear"] == "true" and p["spinorbit"] == "true"
    row = p["atoms"].splitlines()[0].split()
    assert len(row) == 10
    assert [float(x) for x in row[7:]] == [0.2, 90.0, 45.0]
    assert any("fully relativistic" in w for w in notes.warnings)


def test_starting_magnetization_is_halved_and_species_labels_map_to_elements():
    species = ("ATOMIC_SPECIES\n  Si1 28.0 Si.header-only.UPF\n  Si2 28.0 Si.header-only.UPF\n")
    pos = "ATOMIC_POSITIONS crystal\n  Si1 0 0 0\n  Si2 0.25 0.25 0.25\n"
    text = mini(species=species, positions=pos,
                system="  nspin = 2, starting_magnetization(1) = 1.0\n"
                       "  starting_magnetization(2) = -0.6\n"
                       "  occupations = 'smearing', degauss = 0.01\n").replace("ntyp = 1", "ntyp = 2")
    p, _ = run(text)
    rows = [l.split() for l in p["atoms"].splitlines()]
    assert [r[0] for r in rows] == ["Si", "Si"]
    assert [float(r[5]) for r in rows] == [0.5, -0.3]
    assert p["pseudopotential"].split() == ["Si", "Si.header-only.UPF"]


def test_nbnd_without_upf_warns(tmp_path):
    text = mini(system="  nbnd = 12\n").replace("'pseudo'", "'%s'" % str(tmp_path))
    p, notes = run(text)
    assert "unoccupied_states_per_kpoint" not in p
    assert any("nbnd" in w for w in notes.warnings)


def test_charge_and_symmetry_and_vdw():
    p, _ = run(mini(system="  tot_charge = 1.0, nosym = .true., noinv = .true.\n"
                           "  vdw_corr = 'grimme-d3', dftd3_version = 4\n  nbnd = 6\n"))
    assert float(p["system_charge"]) == 1.0
    assert p["use_symmetry"] == "0" and p["time_reversal"] == "false"
    assert p["vdw_corr"] == "DFT-D3" and p["dftd3_version"] == "4"
    assert p["unoccupied_states_per_kpoint"] == "2"     # 7 electrons -> 4 occupied states


# ---------------------------------------------------------------------------------------------
# Calculation types and k-points
# ---------------------------------------------------------------------------------------------


def test_calculation_modes():
    p, _ = run(mini(control="  calculation = 'relax', forc_conv_thr = 2.0d-3, nstep = 30\n"))
    assert p["calculation_mode"] == "Relax Structure"
    assert float(p["relax_max_force"]) == pytest.approx(1e-3)      # Ry/bohr -> Ha/bohr
    assert p["max_md_steps"] == "30" and p["relax_method"] == "LBFGS"
    p, _ = run(mini(control="  calculation = 'vc-relax'\n"))
    assert p["cell_relax"] == "true"
    p, _ = run(mini(control="  calculation = 'nscf'\n"))
    assert p["calculation_mode"] == "NSCF"
    p, notes = run(mini(control="  calculation = 'md'\n"))
    assert p["calculation_mode"] == "Constant Volume And Energy"
    assert any("time step" in w for w in notes.warnings)
    with pytest.raises(ConversionError):
        run(mini(control="  calculation = 'bands'\n"))


def test_kpoint_variants():
    p, _ = run(mini(kpoints="K_POINTS gamma\n"))
    assert p["kpoint_mesh"] == "1 1 1" and p["kpoint_is_shift"] == "0 0 0"
    p, _ = run(mini(kpoints="K_POINTS tpiba\n2\n 0 0 0 1\n 0.5 0.5 0.5 3\n"))
    k = [[float(x) for x in l.split()] for l in p["kpoints"].splitlines()]
    assert np.allclose(k, [[0, 0, 0, 1], [0.5, 0.5, 0.5, 3]])   # simple cubic: tpiba = crystal
    assert p["kpoint_mesh"].split()[0] == "-1"
    with pytest.raises(ConversionError):
        run(mini(kpoints="K_POINTS crystal_b\n2\n 0 0 0 10\n 0.5 0 0 1\n"))
    with pytest.raises(ConversionError):
        run(mini(kpoints="K_POINTS automatic\n 2 2 2 0 0 2\n"))


def test_explicit_kpoints_follow_the_rmg_basis_change():
    """For QE's bcc vectors RMG changes the primitive basis; k . r must be unchanged."""
    pos = "ATOMIC_POSITIONS crystal\n  Si 0 0 0\n  Si 0.1 0.2 0.35\n"
    kp = "K_POINTS crystal\n1\n 0.25 0.1 0.4 1\n"
    p, _ = run(mini(ibrav=3, system="  celldm(1) = 10.0\n", positions=pos, kpoints=kp))
    st = deck_structure(p)
    kf = np.array([float(x) for x in p["kpoints"].split()[:3]])
    phase_in = np.array([0.25, 0.1, 0.4]) @ np.array([0.1, 0.2, 0.35])
    assert kf @ st.frac[1] == pytest.approx(phase_in, abs=1e-10)


# ---------------------------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("extra", [
    dict(system="  lda_plus_u = .true., Hubbard_U(1) = 4.0\n"),
    dict(control="  tefield = .true.\n"),
    dict(system="  assume_isolated = 'mt'\n"),
    dict(system="  occupations = 'from_input'\n"),
    dict(system="  vdw_corr = 'ts'\n"),
    dict(system="  input_dft = 'pbe0'\n", kpoints="K_POINTS crystal\n1\n 0 0 0 1\n"),
])
def test_unsupported_physics_is_refused(extra):
    with pytest.raises(ConversionError):
        run(mini(**extra))


def test_unsupported_cards_are_refused():
    with pytest.raises(ConversionError):
        run(mini() + "HUBBARD {ortho-atomic}\nU Si-3p 1.0\n")
    with pytest.raises(ConversionError):
        run(mini() + "CONSTRAINTS\n1\n'distance' 1 2\n")


def test_upf_v1_is_refused(tmp_path):
    (tmp_path / "Si.v1.UPF").write_text(
        "<PP_HEADER>\n  Si  Element\n  NC  Norm - Conserving\n  4.0  Z valence\n</PP_HEADER>\n")
    text = mini(species="ATOMIC_SPECIES\n  Si 28.0 Si.v1.UPF\n")
    with pytest.raises(ConversionError):
        run(text, pseudo_dir=str(tmp_path))
    p, _ = run(text, pseudo_dir=str(tmp_path), lenient=True)
    assert p["pseudo_dir"] == str(tmp_path)


def test_unknown_keys_warn():
    p, notes = run(mini(system="  some_new_flag = 3\n"))
    assert any("some_new_flag" in w for w in notes.warnings)


# ---------------------------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------------------------


def test_cli_writes_deck_and_notes(tmp_path, capsys):
    out = tmp_path / "rmg.in"
    rc = main([os.path.join(DATA, "zgnr_hse.pw.in"), "-o", str(out)])
    assert rc == 0
    deck = parse_rmg_deck(out.read_text())
    assert deck["exchange_correlation_type"] == "hse"
    notes = (tmp_path / "rmg.in.notes").read_text()
    assert "exx_qpoint_mesh" in notes and "NOTE:" in notes
    assert "qe2rmg: note:" in capsys.readouterr().err


def test_cli_grid_override_and_refusal(tmp_path, capsys):
    out = tmp_path / "si.in"
    assert main([os.path.join(DATA, "si_diamond.pw.in"), "-o", str(out),
                 "--grid", "20", "20", "20"]) == 0
    assert parse_rmg_deck(out.read_text())["wavefunction_grid"] == "20 20 20"
    assert main([os.path.join(DATA, "si_diamond.pw.in"), "-o", str(out),
                 "--grid-spacing", "0.3"]) == 0
    assert parse_rmg_deck(out.read_text())["wavefunction_grid"] == "28 28 28"
    bad = tmp_path / "bad.in"
    bad.write_text(mini(control="  calculation = 'bands'\n"))
    assert main([str(bad), "-o", str(tmp_path / "x.in"), "--pseudo-dir", PSEUDO]) == 2
    assert "error" in capsys.readouterr().err
