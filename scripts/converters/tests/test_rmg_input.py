"""Tests for the shared RMG deck model, grids, lattice orientation and UPF reader."""

import math
import os

import numpy as np
import pytest

from conftest import PSEUDO

from rmg_input import (BOHR_PER_ANGSTROM, ConversionError, Notes, RmgDeck, Structure,
                       check_upf_for_rmg, clamp_spin_ratio, deck_structure, element_from_label,
                       equivalent_cutoff_ry, min_points_for_cutoff, next_grid_size,
                       orient_for_rmg, parse_rmg_deck, planar_anisotropy,
                       potential_grid_refinement, read_upf_header, rmg_lat2ibrav, rmg_rotate,
                       wavefunction_grid, write_structure)


# ---------------------------------------------------------------------------------------------
# Deck model
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("bad", ["a=b", 'say "hi"', "x # comment"])
def test_deck_rejects_values_rmg_cannot_parse(bad):
    d = RmgDeck()
    with pytest.raises(ValueError):
        d.set("description", bad)


def test_deck_rejects_bad_keys_and_comment_lines_in_blocks():
    d = RmgDeck()
    with pytest.raises(ValueError):
        d.set("bad key", "1")
    with pytest.raises(ValueError):
        d.set("1abc", "1")
    with pytest.raises(ValueError):
        d.set_block("atoms", ["  Si 0 0 0 1", "# not allowed"])
    with pytest.raises(ValueError):
        d.set("description", "two\nlines")


def test_deck_render_parses_back_exactly():
    d = RmgDeck()
    d.set("description", "Si test")
    d.set("spin_polarization", True)
    d.set("max_scf_steps", 40)
    d.set("charge_density_mixing", 0.3)
    d.set("kpoint_mesh", [4, 4, 4])
    d.set_block("atoms", ["  Si 0.0 0.0 0.0 1", "  Si 0.25 0.25 0.25 1"])
    text = d.render()
    assert not any(l.lstrip().startswith("#") for l in text.splitlines())
    p = parse_rmg_deck(text)
    assert p["spin_polarization"] == "true"
    assert p["max_scf_steps"] == "40"
    assert float(p["charge_density_mixing"]) == 0.3
    assert p["kpoint_mesh"] == "4 4 4"
    assert p["atoms"].splitlines()[1].split() == ["Si", "0.25", "0.25", "0.25", "1"]
    assert list(p) == d.keys()


def test_parse_rmg_deck_follows_rmg_comment_rules():
    text = '# a comment\nkpoint_mesh = "2 2 2"  # trailing\natoms = "\n  H 0 0 0 1\n"\n'
    p = parse_rmg_deck(text)
    assert p == {"kpoint_mesh": "2 2 2", "atoms": "H 0 0 0 1"}
    with pytest.raises(ValueError):
        parse_rmg_deck('kpoint_mesh "2 2 2"')


def test_write_structure_columns():
    lat = np.diag([10.0, 10.0, 10.0])
    st = Structure(lat, ["Fe", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]],
                   movable=[[1, 1, 1], [1, 0, 1]], spin=[(0.3,), (-0.1,)])
    d = RmgDeck()
    write_structure(d, st, spin_columns=1)
    p = parse_rmg_deck(d.render())
    rows = [l.split() for l in p["atoms"].splitlines()]
    assert rows[0][4:] == ["1", "1", "1", "+0.300000"]
    assert rows[1][4:] == ["1", "0", "1", "-0.100000"]
    back = deck_structure(p)
    assert np.allclose(back.lattice_bohr, lat, atol=1e-10)
    assert back.movable.tolist() == [[1, 1, 1], [1, 0, 1]]
    assert back.moments == [(0.3,), (-0.1,)]
    d2 = RmgDeck()
    write_structure(d2, Structure(lat, ["Fe"], [[0, 0, 0]], spin=[(0.2, 45.0, 90.0)]), 3)
    row = parse_rmg_deck(d2.render())["atoms"].split()
    assert len(row) == 10 and float(row[8]) == 45.0 and float(row[9]) == 90.0


def test_clamp_spin_ratio():
    n = Notes()
    assert clamp_spin_ratio(0.2, n) == 0.2
    assert clamp_spin_ratio(0.8, n) == 0.5
    assert clamp_spin_ratio(-0.9, n) == -0.5
    assert len(n.warnings) == 2


def test_notes_unsupported_refuses_unless_lenient():
    with pytest.raises(ConversionError):
        Notes().unsupported("x")
    n = Notes(lenient=True)
    n.unsupported("x")
    assert n.warnings and n.warnings[0].startswith("DROPPED")


# ---------------------------------------------------------------------------------------------
# Grids
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("n,expected", [(1, 4), (13, 16), (17, 20), (21, 24), (28, 28),
                                        (41, 48), (89, 96), (97, 100)])
def test_next_grid_size(n, expected):
    assert next_grid_size(n) == expected


def test_equivalent_cutoff_matches_rmg_for_a_known_deck():
    # The 10-ZGNR HSE06 production deck: 768 352 28 on 67.2 x 31 x 2.4663 A is reported by
    # RMG as 356.2 Ry (the smallest spacing direction limits the sphere).
    lat = np.diag([67.2, 31.0, 2.4663]) * BOHR_PER_ANGSTROM
    assert equivalent_cutoff_ry(lat, (768, 352, 28)) == pytest.approx(356.2, abs=0.1)
    assert min_points_for_cutoff(2.4663 * BOHR_PER_ANGSTROM, 356.2) == 28


def test_grid_from_cutoff_is_sufficient_and_minimal():
    lat = np.array([[0.0, 5.13, 5.13], [5.13, 0.0, 5.13], [5.13, 5.13, 0.0]])
    for ecut in (20.0, 30.0, 45.0, 80.0, 120.0):
        g = wavefunction_grid(lat, ecut_ry=ecut)
        assert g.ecut_equiv_ry >= ecut - 1e-9
        n = g.grid[0]
        assert all(x == n for x in g.grid) and n % 4 == 0
        smaller = [m for m in range(4, n, 4) if next_grid_size(m) == m]
        if smaller:
            assert equivalent_cutoff_ry(lat, (smaller[-1],) * 3) < ecut


def test_grid_spacing_and_explicit_grid():
    lat = np.diag([10.0, 12.0, 20.0])
    g = wavefunction_grid(lat, spacing_bohr=0.25)
    assert g.grid == (40, 48, 80)
    n = Notes()
    g2 = wavefunction_grid(lat, grid=(30, 36, 60), notes=n)
    assert g2.grid == (30, 36, 60) and not n.warnings
    with pytest.raises(ConversionError):
        wavefunction_grid(lat)


def test_grid_respects_rmg_planar_anisotropy_limit():
    lat = np.diag([2.0, 2.0, 40.0])
    n = Notes()
    g = wavefunction_grid(lat, spacing_bohr=0.45, notes=n)
    assert planar_anisotropy(lat, g.grid) <= 1.5
    assert any("anisotropy" in w for w in n.warnings)
    n2 = Notes()
    wavefunction_grid(lat, grid=(8, 8, 96), notes=n2)
    assert any("anisotropy" in w for w in n2.warnings)


@pytest.mark.parametrize("dual,ratio", [(4, 2), (2, 2), (8, 3), (9, 3), (12, 4), (16, 4), (25, 4)])
def test_potential_grid_refinement(dual, ratio):
    assert potential_grid_refinement(dual * 30.0, 30.0) == ratio


# ---------------------------------------------------------------------------------------------
# Lattice orientation (RMG Lattice::lat2ibrav / rotate_vectors)
# ---------------------------------------------------------------------------------------------


def _distances(lat, frac, rmax=4.0):
    """Sorted interatomic distances (periodic images included) up to ``rmax`` bohr."""
    cart = np.asarray(frac) @ lat
    shifts = np.array(np.meshgrid(*[range(-3, 4)] * 3)).T.reshape(-1, 3) @ lat
    out = []
    for ci in cart:
        for cj in cart:
            d = np.linalg.norm(cj + shifts - ci, axis=1)
            out.extend(d[d < rmax])
    return np.sort(np.round(out, 9))


def test_fcc_is_rotated_with_unchanged_fractions():
    a = 10.26
    qe = np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]]) * a / 2
    assert rmg_lat2ibrav(qe) == 2
    f = np.array([[0, 0, 0], [0.25, 0.25, 0.25]])
    o = orient_for_rmg(qe, f)
    assert o.changed and not o.mirrored
    assert np.allclose(o.lattice, np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]]) * a / 2)
    assert np.allclose(o.frac, f)
    assert np.allclose(o.lattice @ o.lattice.T, qe @ qe.T)


def test_bcc_qe_vectors_are_rebased_and_the_crystal_is_kept():
    a = 5.42
    qe = np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]]) * a / 2
    f = np.array([[0.0, 0.0, 0.0], [0.1, 0.2, 0.35]])
    o = orient_for_rmg(qe, f)
    assert o.changed and not np.allclose(o.basis_change, np.eye(3))
    # the written cell is a fixed point of RMG's re-orientation
    assert np.allclose(rmg_rotate(o.lattice), o.lattice)
    assert abs(np.linalg.det(o.lattice)) == pytest.approx(abs(np.linalg.det(qe)))
    assert np.allclose(_distances(o.lattice, o.frac), _distances(qe, f))


def test_rotated_orthorhombic_and_hexagonal_cells():
    ortho = np.array([[0.0, 7.0, 0.0], [-5.0, 0.0, 0.0], [0.0, 0.0, 9.0]])
    f = np.array([[0.1, 0.2, 0.3], [0.6, 0.4, 0.9]])
    o = orient_for_rmg(ortho, f)
    assert np.allclose(o.lattice, np.diag([7.0, 5.0, 9.0]))
    assert np.allclose(_distances(o.lattice, o.frac), _distances(ortho, f))
    hexa = np.array([[1, 0, 0], [-0.5, math.sqrt(3) / 2, 0], [0, 0, 1.6]]) * 4.6
    o2 = orient_for_rmg(hexa, f)
    assert not o2.changed and np.allclose(o2.frac, f)


def test_body_centred_tetragonal_is_refused():
    a, c = 5.0, 7.0
    bct = np.array([[a / 2, -a / 2, c / 2], [a / 2, a / 2, c / 2], [-a / 2, -a / 2, c / 2]])
    with pytest.raises(ConversionError):
        orient_for_rmg(bct, [[0, 0, 0]])


def test_triclinic_is_left_alone():
    tri = np.array([[5.0, 0.1, 0.2], [0.3, 6.0, 0.1], [0.2, 0.4, 7.0]])
    o = orient_for_rmg(tri, [[0.1, 0.2, 0.3]])
    assert not o.changed and np.allclose(o.lattice, tri)


# ---------------------------------------------------------------------------------------------
# Pseudopotentials and labels
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("label,el", [("Fe1", "Fe"), ("Fe_up", "Fe"), ("Fe_pv", "Fe"),
                                      ("C2", "C"), ("Co", "Co"), ("C_h", "C"), ("H", "H"),
                                      ("Si", "Si"), ("O1", "O"), ("FE", "Fe"), ("Xx", None)])
def test_element_from_label(label, el):
    assert element_from_label(label) == el


def test_read_upf_v2_header():
    h = read_upf_header(os.path.join(PSEUDO, "Fe.header-only.UPF"))
    assert (h.version, h.element, h.z_valence, h.pseudo_type, h.has_so) == (2, "Fe", 16.0, "NC",
                                                                           False)
    check_upf_for_rmg(h, Notes())


def test_upf_v1_and_paw_are_reported(tmp_path):
    v1 = tmp_path / "Si.v1.UPF"
    v1.write_text("<PP_INFO>\n</PP_INFO>\n<PP_HEADER>\n   0   Version Number\n  Si   Element\n"
                  "   NC  Norm - Conserving pseudopotential\n    F  Nonlinear Core Correction\n"
                  " SLA PW PBE PBE  PBE  Exchange-Correlation functional\n"
                  "    4.00000000000  Z valence\n</PP_HEADER>\n")
    h = read_upf_header(str(v1))
    assert (h.version, h.element, h.z_valence, h.pseudo_type) == (1, "Si", 4.0, "NC")
    with pytest.raises(ConversionError):
        check_upf_for_rmg(h, Notes())
    paw = tmp_path / "Si.paw.UPF"
    paw.write_text('<UPF version="2.0.1"><PP_HEADER element="Si" pseudo_type="PAW" '
                   'z_valence="4.0"/></UPF>')
    with pytest.raises(ConversionError):
        check_upf_for_rmg(read_upf_header(str(paw)), Notes())
