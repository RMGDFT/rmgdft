"""Tests for the VASP (POSCAR/INCAR/KPOINTS) -> RMG converter."""

import math
import os

import numpy as np
import pytest

from conftest import DATA, PSEUDO, read_data
from rmg_input import BOHR_PER_ANGSTROM, HA_EV, ConversionError, deck_structure, parse_rmg_deck
from vasp2rmg import (convert, expand_magmom, main, mesh_from_kspacing, mesh_from_length,
                      parse_incar, parse_kpoints, parse_poscar, parse_potcar)

SI_POSCAR = read_data("vasp_si", "POSCAR")
SI_INCAR = read_data("vasp_si", "INCAR")
SI_KPOINTS = read_data("vasp_si", "KPOINTS")
FE_POSCAR = read_data("vasp_fe", "POSCAR")
FE_INCAR = read_data("vasp_fe", "INCAR")
FE_KPOINTS = read_data("vasp_fe", "KPOINTS")
FE_UPF = os.path.join(PSEUDO, "Fe.header-only.UPF")
SI_UPF = os.path.join(PSEUDO, "Si.header-only.UPF")
GAMMA4 = "mesh\n0\nGamma\n4 4 4\n"


def run(poscar=SI_POSCAR, incar=SI_INCAR, kpoints=SI_KPOINTS, **kw):
    deck, notes = convert(poscar, incar, kpoints, **kw)
    return parse_rmg_deck(deck.render()), notes


def fake_potcar(*species):
    """POTCAR-shaped text with only the header lines the converter reads (no PAW data)."""
    out = []
    for title, zval, enmax in species:
        out.append("  PAW_PBE %s 01Jan2000\n   TITEL  = PAW_PBE %s 01Jan2000\n"
                   "   POMASS =   1.0; ZVAL   =   %.3f    mass and valenz\n"
                   "   ENMAX  =  %.3f; ENMIN  =  100.000 eV\n End of Dataset\n"
                   % (title, title, zval, enmax))
    return "".join(out)


# ---------------------------------------------------------------------------------------------
# POSCAR
# ---------------------------------------------------------------------------------------------


def test_poscar_vasp5_direct():
    p = parse_poscar(SI_POSCAR)
    assert p.symbols == ["Si", "Si"]
    assert np.allclose(p.lattice, np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]) * 5.43)
    assert np.allclose(p.frac, [[0, 0, 0], [0.25, 0.25, 0.25]])
    assert p.selective is None


def test_poscar_vasp4_needs_species():
    v4 = "Si\n5.43\n0 .5 .5\n.5 0 .5\n.5 .5 0\n2\nDirect\n0 0 0\n.25 .25 .25\n"
    with pytest.raises(ConversionError):
        parse_poscar(v4)
    assert parse_poscar(v4, ["Si"]).symbols == ["Si", "Si"]
    with pytest.raises(ConversionError):
        parse_poscar(SI_POSCAR, ["Ge"])


def test_poscar_cartesian_negative_scale_and_three_scales():
    cart = ("x\n-24.0\n1 0 0\n0 2 0\n0 0 3\nO\n1\nCartesian\n0.5 0.5 0.5\n")
    p = parse_poscar(cart)
    assert abs(np.linalg.det(p.lattice)) == pytest.approx(24.0)        # negative = volume
    s = (24.0 / 6.0) ** (1 / 3)
    assert np.allclose(p.frac @ p.lattice, [[0.5 * s, 0.5 * s, 0.5 * s]])
    three = "x\n1.0 2.0 3.0\n1 0 0\n0 1 0\n0 0 1\nO\n1\nCartesian\n0.5 0.5 0.5\n"
    p3 = parse_poscar(three)
    assert np.allclose(p3.lattice, np.diag([1.0, 2.0, 3.0]))
    assert np.allclose(p3.frac, [[0.5, 0.5, 0.5]])


def test_poscar_selective_dynamics_and_species_suffixes():
    txt = ("x\n1.0\n4 0 0\n0 4 0\n0 0 4\nFe_pv O_s\n1 1\nSelective dynamics\nDirect\n"
           "0 0 0 F F F\n0.5 0.5 0.5 T F T\n")
    p = parse_poscar(txt)
    assert p.symbols == ["Fe", "O"]
    assert p.selective.tolist() == [[0, 0, 0], [1, 0, 1]]
    d, _ = run(poscar=txt, incar="ENCUT = 400\nIBRION = 2\nNSW = 10\n", kpoints=GAMMA4)
    rows = [l.split() for l in d["atoms"].splitlines()]
    assert rows[0][4:7] == ["0", "0", "0"] and rows[1][4:7] == ["1", "0", "1"]


# ---------------------------------------------------------------------------------------------
# INCAR, MAGMOM, KPOINTS, POTCAR
# ---------------------------------------------------------------------------------------------


def test_incar_parsing():
    t = parse_incar("ENCUT = 520 # cutoff\nISMEAR = 0; SIGMA = 0.05 ! smearing\n"
                    "MAGMOM = 2*3.0 \\\n  1*-1\nSYSTEM = \"my run\"\n")
    assert t == {"ENCUT": "520", "ISMEAR": "0", "SIGMA": "0.05", "MAGMOM": "2*3.0    1*-1",
                 "SYSTEM": "my run"}


@pytest.mark.parametrize("txt,val", [("2*3.0 1*-1 0", [3.0, 3.0, -1.0, 0.0]),
                                     ("4*0", [0.0] * 4), ("1 -1", [1.0, -1.0])])
def test_magmom_expansion(txt, val):
    assert expand_magmom(txt) == val


def test_kpoint_shifts():
    lat = np.eye(3) * 4.0
    assert parse_kpoints("a\n0\nGamma\n4 4 4\n", lat)[1:3] == ((4, 4, 4), (0, 0, 0))
    assert parse_kpoints("a\n0\nMonkhorst-Pack\n4 3 2\n", lat)[2] == (1, 0, 1)
    assert parse_kpoints("a\n0\nMonkhorst\n3 3 3\n0.5 0 0\n", lat)[2] == (1, 0, 0)
    assert parse_kpoints("a\n0\nGamma\n2 2 2\n0.5 0.5 0\n", lat)[2] == (1, 1, 0)
    assert parse_kpoints("a\n0\nM\n4 4 4\n0.5 0.5 0.5\n", lat)[2] == (0, 0, 0)
    with pytest.raises(ConversionError):
        parse_kpoints("a\n0\nGamma\n4 4 4\n0.25 0 0\n", lat)
    with pytest.raises(ConversionError):
        parse_kpoints("a\n10\nLine-mode\nReciprocal\n0 0 0 G\n0.5 0 0 X\n", lat)


def test_fully_automatic_and_kspacing_meshes():
    lat = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]) * 5.43
    b = np.linalg.norm(np.linalg.inv(lat).T, axis=1)          # 1/A, no 2 pi
    assert mesh_from_length(lat, 20.0) == tuple(max(1, int(20 * x + 0.5)) for x in b)
    kp = parse_kpoints("auto\n0\nAuto\n20\n", lat)
    assert kp[0] == "mesh" and kp[2] == (0, 0, 0)
    assert mesh_from_kspacing(lat, 0.3) == tuple(math.ceil(2 * math.pi * x / 0.3) for x in b)
    d, notes = run(kpoints=None, incar=SI_INCAR + "KSPACING = 0.3\nKGAMMA = .FALSE.\n")
    m = [int(x) for x in d["kpoint_mesh"].split()]
    assert m == list(mesh_from_kspacing(lat, 0.3))
    assert d["kpoint_is_shift"] == " ".join("1" if x % 2 == 0 else "0" for x in m)


def test_explicit_reciprocal_kpoints():
    d, _ = run(kpoints="list\n2\nReciprocal\n0 0 0 1\n0.5 0.5 0 2\n")
    k = [[float(x) for x in l.split()] for l in d["kpoints"].splitlines()]
    assert np.allclose(k, [[0, 0, 0, 1], [0.5, 0.5, 0, 2]])     # fcc keeps its basis (M = I)


def test_potcar_reader():
    pc = parse_potcar(fake_potcar(("Fe_pv", 14.0, 293.2), ("O", 6.0, 400.0)))
    assert pc.elements == ["Fe", "O"] and pc.zval == [14.0, 6.0] and pc.enmax == [293.2, 400.0]


# ---------------------------------------------------------------------------------------------
# Fixtures: round trips
# ---------------------------------------------------------------------------------------------


def test_si_fixture():
    d, notes = run()
    st = deck_structure(d)
    lat = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]) * 5.43 * BOHR_PER_ANGSTROM
    assert np.allclose(st.lattice_bohr @ st.lattice_bohr.T, lat @ lat.T, atol=1e-8)
    assert np.allclose(st.frac, [[0, 0, 0], [0.25, 0.25, 0.25]], atol=1e-10)
    assert d["occupations_type"] == "MethfesselPaxton" and d["MP_order"] == "0"
    assert float(d["occupation_electron_temperature_eV"]) == 0.05
    assert d["kpoint_mesh"] == "4 4 4" and d["kpoint_is_shift"] == "1 1 1"
    assert float(d["energy_convergence_criterion"]) == pytest.approx(1e-6 / HA_EV)
    assert float(d["rms_convergence_criterion"]) == pytest.approx(1e-3 * (2e-6 / HA_EV) ** 0.5)
    assert d["max_scf_steps"] == "80"
    g = [int(x) for x in d["wavefunction_grid"].split()]
    h = np.linalg.norm(st.lattice_bohr, axis=1) / g
    assert np.all(h <= math.pi / math.sqrt(300.0 / (HA_EV / 2)) + 1e-12)
    assert "pseudopotential" not in d
    assert any("PAW" in w for w in notes.warnings)
    assert any("NBANDS" in w for w in notes.warnings)
    d2, _ = run(pseudos={"Si": SI_UPF})
    assert d2["unoccupied_states_per_kpoint"] == "4"
    assert d2["pseudopotential"].split() == ["Si", "Si.header-only.UPF"]
    assert d2["pseudo_dir"] == PSEUDO


def test_fe_fixture():
    d, notes = run(FE_POSCAR, FE_INCAR, FE_KPOINTS, pseudos={"Fe": FE_UPF})
    st = deck_structure(d)
    assert np.allclose(st.lattice_bohr, np.eye(3) * 2.86 * BOHR_PER_ANGSTROM, atol=1e-8)
    assert np.allclose(st.cart_bohr[1], [1.43 * BOHR_PER_ANGSTROM] * 3, atol=1e-8)
    assert d["spin_polarization"] == "true"
    assert [m[0] for m in st.moments] == [pytest.approx(2.5 / 32, abs=1e-6)] * 2
    assert d["exchange_correlation_type"] == "pbe"
    assert d["occupations_type"] == "MethfesselPaxton" and d["MP_order"] == "1"
    assert float(d["occupation_electron_temperature_eV"]) == 0.1
    assert d["kpoint_mesh"] == "8 8 8" and d["kpoint_is_shift"] == "0 0 0"
    assert float(d["charge_density_mixing"]) == 0.2
    assert d["max_scf_steps"] == "100"
    # ISPIN 2: RMG's energy test includes the magnetization residual, so EDIFF maps onto it directly
    assert float(d["energy_convergence_criterion"]) == pytest.approx(1e-7 / HA_EV)
    assert float(d["rms_convergence_criterion"]) == pytest.approx(1e-3 * (2e-7 / HA_EV) ** 0.5)
    # without valence charges the seed falls back to MAGMOM / 10
    d2, notes2 = run(FE_POSCAR, FE_INCAR, FE_KPOINTS)
    assert [float(l.split()[5]) for l in d2["atoms"].splitlines()] == [0.25, 0.25]
    assert any("MAGMOM / 10" in w for w in notes2.warnings)


def test_magmom_count_and_nupdown():
    with pytest.raises(ConversionError):
        run(FE_POSCAR, FE_INCAR.replace("2*2.5", "2.5"), FE_KPOINTS)
    d, _ = run(FE_POSCAR, FE_INCAR + "NUPDOWN = 4.4\n", FE_KPOINTS)
    assert float(d["tot_magnetization"]) == 4.4
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + "NUPDOWN = 2\n")


def test_noncollinear_magmom_angles():
    d, _ = run(FE_POSCAR, FE_INCAR.replace("ISPIN = 2", "LNONCOLLINEAR = .TRUE.").replace(
        "MAGMOM = 2*2.5", "MAGMOM = 0 0 2  2 0 0"), FE_KPOINTS, pseudos={"Fe": FE_UPF})
    rows = [[float(x) for x in l.split()[1:]] for l in d["atoms"].splitlines()]
    assert d["noncollinear"] == "true"
    assert rows[0][6:] == [pytest.approx(2 / 32, abs=1e-6), 0.0, 0.0]
    assert rows[1][6:] == [pytest.approx(2 / 32, abs=1e-6), 90.0, 0.0]


# ---------------------------------------------------------------------------------------------
# Functionals and exact exchange
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("tags,xc", [("GGA = PE", "pbe"), ("GGA = PS", "pbesol"),
                                     ("GGA = 91", "pw91"), ("GGA = CA", "pz"),
                                     ("METAGGA = SCAN", "scan"), ("METAGGA = R2SCAN", "r2scan"),
                                     ("LHFCALC = .TRUE.", "pbe0"),
                                     ("LHFCALC = .TRUE.\nHFSCREEN = 0.2", "hse"),
                                     ("LHFCALC = .TRUE.\nGGA = B3\nAEXX = 0.2", "b3lyp")])
def test_functional_mapping(tags, xc):
    d, _ = run(incar=SI_INCAR + tags + "\n", kpoints=GAMMA4)
    assert d["exchange_correlation_type"] == xc
    assert ("exx_qpoint_mesh" in d) == ("LHFCALC" in tags)


def test_hybrid_qmesh_from_nkred_and_refusals():
    d, _ = run(incar=SI_INCAR + "LHFCALC = .TRUE.\nHFSCREEN = 0.2\nNKRED = 2\n", kpoints=GAMMA4)
    assert d["exx_qpoint_mesh"] == "2 2 2"
    d, _ = run(incar=SI_INCAR + "LHFCALC = .TRUE.\nAEXX = 0.3\n", kpoints=GAMMA4)
    assert float(d["exx_fraction"]) == 0.3
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + "LHFCALC = .TRUE.\nNKRED = 3\n", kpoints=GAMMA4)
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + "LHFCALC = .TRUE.\nHFSCREEN = 0.3\n", kpoints=GAMMA4)
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + "LHFCALC = .TRUE.\nAGGAX = 0.5\n", kpoints=GAMMA4)


# ---------------------------------------------------------------------------------------------
# Occupations, charge, calculation types
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("ismear,otype,extra", [(-5, "Tetrahedron", ("tetra_method", "Bloechl")),
                                                (-4, "Tetrahedron", ("tetra_method", "Linear")),
                                                (-1, "Fermi Dirac", None),
                                                (2, "MethfesselPaxton", ("MP_order", "2"))])
def test_ismear(ismear, otype, extra):
    d, _ = run(incar="ENCUT = 300\nISMEAR = %d\nSIGMA = 0.1\n" % ismear)
    assert d["occupations_type"] == otype
    if extra:
        assert d[extra[0]] == extra[1]


def test_ismear_minus2_refused():
    with pytest.raises(ConversionError):
        run(incar="ENCUT = 300\nISMEAR = -2\n")


def test_vasp_defaults_without_incar_tags():
    d, _ = run(incar="ENCUT = 300\n")
    assert d["occupations_type"] == "MethfesselPaxton" and d["MP_order"] == "1"
    assert float(d["occupation_electron_temperature_eV"]) == 0.2
    assert d["max_scf_steps"] == "60"


def test_nelect_needs_potcar_or_charge():
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + "NELECT = 7\n")
    pc = fake_potcar(("Si", 4.0, 245.3))
    d, _ = run(incar=SI_INCAR + "NELECT = 7\n", potcar=pc)
    assert float(d["system_charge"]) == 1.0
    d, _ = run(incar=SI_INCAR + "NELECT = 7\n", system_charge=-1.0)
    assert float(d["system_charge"]) == -1.0


def test_potcar_gives_default_encut_and_checks_species():
    pc = fake_potcar(("Si", 4.0, 245.3))
    d, notes = run(incar="ISMEAR = 0\n", potcar=pc)
    assert any("245.3" in i for i in notes.infos)
    with pytest.raises(ConversionError):
        run(incar="ISMEAR = 0\n")
    with pytest.raises(ConversionError):
        run(potcar=fake_potcar(("Ge", 4.0, 174.0)))


def test_calculation_modes():
    d, _ = run(incar=SI_INCAR + "IBRION = 2\nNSW = 40\nEDIFFG = -0.01\n")
    assert d["calculation_mode"] == "Relax Structure" and d["max_md_steps"] == "40"
    assert float(d["relax_max_force"]) == pytest.approx(0.01 / HA_EV / BOHR_PER_ANGSTROM)
    d, _ = run(incar=SI_INCAR + "IBRION = 1\nNSW = 5\nISIF = 3\n")
    assert d["cell_relax"] == "true"
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + "IBRION = 2\nNSW = 5\nISIF = 4\n")
    d, notes = run(incar=SI_INCAR + "IBRION = 0\nNSW = 100\nPOTIM = 1.0\n")
    assert d["calculation_mode"] == "Constant Volume And Energy"
    d, _ = run(incar=SI_INCAR + "IBRION = 0\nNSW = 100\nSMASS = 0\nTEBEG = 500\n")
    assert d["calculation_mode"] == "Constant Temperature And Energy"
    assert float(d["md_temperature"]) == 500.0
    d, _ = run(incar=SI_INCAR + "ICHARG = 11\n")
    assert d["calculation_mode"] == "NSCF"
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + "IBRION = 6\nNSW = 1\n")


@pytest.mark.parametrize("tag", ["LDAU = .TRUE.", "LEPSILON = .TRUE.", "LDIPOL = .TRUE.",
                                 "IVDW = 4", "SAXIS = 1 0 0\nLSORBIT = .TRUE."])
def test_unsupported_tags_refused(tag):
    with pytest.raises(ConversionError):
        run(incar=SI_INCAR + tag + "\n")


def test_vdw_symmetry_and_unknown_tags():
    d, notes = run(incar=SI_INCAR + "IVDW = 12\nISYM = -1\nFOOBAR = 3\n")
    assert d["vdw_corr"] == "DFT-D3" and d["dftd3_version"] == "4"
    assert d["use_symmetry"] == "0" and d["time_reversal"] == "false"
    assert any("FOOBAR" in w for w in notes.warnings)


def test_pseudo_mapping_must_cover_every_element():
    two = "x\n1.0\n4 0 0\n0 4 0\n0 0 4\nFe Si\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n"
    with pytest.raises(ConversionError):
        run(poscar=two, incar="ENCUT = 300\n", kpoints=GAMMA4, pseudos={"Fe": FE_UPF})
    with pytest.raises(ConversionError):
        run(poscar=two, incar="ENCUT = 300\n", kpoints=GAMMA4,
            pseudos={"Fe": SI_UPF, "Si": SI_UPF})


# ---------------------------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------------------------


def test_cli(tmp_path, capsys):
    out = tmp_path / "fe.rmg"
    rc = main(["-d", os.path.join(DATA, "vasp_fe"), "-o", str(out), "--pseudo", "Fe=" + FE_UPF])
    assert rc == 0
    d = parse_rmg_deck(out.read_text())
    assert d["spin_polarization"] == "true"
    assert "WARNING" in (tmp_path / "fe.rmg.notes").read_text()
    assert "vasp2rmg: warning:" in capsys.readouterr().err
    assert main(["-d", str(tmp_path)]) == 2       # no POSCAR there
