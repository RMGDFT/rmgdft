#!/usr/bin/env python3
"""Convert a Quantum ESPRESSO pw.x input file into an RMG (rmg-dft) input file.

Usage::

    qe2rmg.py pw.in -o rmg.in [--pseudo-dir DIR] [--grid NX NY NZ | --grid-spacing H]
                               [--grid-multiple N] [--notes FILE] [--lenient]

The deck goes to ``-o`` (or stdout). Everything the converter changed, dropped or could not
check is printed on stderr and written to ``<output>.notes``; RMG decks cannot carry comment
lines. Features that would change the physics and that RMG cannot represent make the
converter stop with exit status 2, unless ``--lenient`` is given, in which case they are
dropped and listed as warnings.

The wavefunction cutoff maps onto RMG's real-space grid: a grid of N points along a lattice
vector of length L resolves plane waves up to a kinetic energy of (pi N / L)^2 Ry, so each
dimension gets N >= L sqrt(ecutwfc[Ry]) / pi, rounded up to a multiple of 4 with only
2, 3, 5, 7 as prime factors. The density cutoff maps onto ``potential_grid_refinement``.
See README.md for the complete keyword table.

Library use::

    from qe2rmg import convert
    deck, notes = convert(open("pw.in").read(), input_dir=".")
    print(deck.render())
"""

from __future__ import annotations

import argparse
import ast
import math
import operator
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rmg_input import (BOHR_PER_ANGSTROM, RY_EV, ConversionError, Notes, RmgDeck,  # noqa: E402
                       Structure, UpfHeader, apply_orientation, check_exx_qmesh,
                       check_upf_for_rmg, clamp_spin_ratio, element_from_label, nocc_rmg,
                       potential_grid_refinement, read_upf_header, sanitize_description,
                       set_scf_thresholds, wavefunction_grid, write_kpoint_list, write_outputs,
                       write_structure)

PROG = "qe2rmg"

# ---------------------------------------------------------------------------------------------
# pw.x input parsing
# ---------------------------------------------------------------------------------------------

CARD_NAMES = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS",
              "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES", "ATOMIC_VELOCITIES",
              "ADDITIONAL_K_POINTS", "SOLVENTS", "HUBBARD", "TOTAL_CHARGE")
_CARD_RE = re.compile(r"^\s*(%s)\b\s*[\{\(]?\s*([A-Za-z_\-]*)\s*[\}\)]?\s*$" % "|".join(CARD_NAMES),
                      re.I)


@dataclass
class Card:
    name: str
    option: str
    lines: List[str] = field(default_factory=list)


@dataclass
class PwInput:
    namelists: Dict[str, Dict[str, Any]]
    cards: Dict[str, Card]

    def nl(self, name: str) -> Dict[str, Any]:
        return self.namelists.get(name, {})


def _strip_comment(line: str) -> str:
    """Remove a '!' comment that is not inside quotes."""
    q = None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "!":
            return line[:i]
    return line


def fortran_value(tok: str) -> Any:
    """Convert one Fortran namelist token to a Python value."""
    t = tok.strip()
    if len(t) >= 2 and t[0] == t[-1] and t[0] in "'\"":
        return t[1:-1].replace(t[0] * 2, t[0])
    low = t.lower()
    if re.fullmatch(r"\.?(t|true)\.?", low):
        return True
    if re.fullmatch(r"\.?(f|false)\.?", low):
        return False
    if re.fullmatch(r"[+-]?\d+", t):
        return int(t)
    try:
        return float(re.sub(r"[dD]", "e", t))
    except ValueError:
        return t


_TOKEN_RE = re.compile(r"""'(?:[^']|'')*'|"(?:[^"]|"")*"|=|,|[^\s,=()'"]+(?:\([^)]*\))?""")


def _parse_namelist_body(body: str) -> Dict[str, Any]:
    toks = _TOKEN_RE.findall(body)
    out: Dict[str, Any] = {}
    i = 0
    while i < len(toks):
        if i + 1 < len(toks) and toks[i + 1] == "=":
            key = re.sub(r"\s+", "", toks[i]).lower()
            j = i + 2
            vals = []
            while j < len(toks) and not (j + 1 < len(toks) and toks[j + 1] == "="):
                if toks[j] != ",":
                    vals.append(fortran_value(toks[j]))
                j += 1
            if not vals:
                raise ConversionError("namelist entry %r has no value" % key)
            m = re.fullmatch(r"([a-z_0-9]+)\((\d+)\)", key)
            if len(vals) == 1:
                out[key] = vals[0]
            else:
                base, start = (m.group(1), int(m.group(2))) if m else (key, 1)
                for k, v in enumerate(vals):
                    out["%s(%d)" % (base, start + k)] = v
            i = j
        elif toks[i] == ",":
            i += 1
        else:
            raise ConversionError("cannot parse namelist text near %r" % toks[i])
    return out


def parse_pw_input(text: str) -> PwInput:
    """Parse namelists (&CONTROL, &SYSTEM, ...) and cards of a pw.x input."""
    lines = [_strip_comment(l) for l in text.splitlines()]
    joined = "\n".join(lines)
    namelists: Dict[str, Dict[str, Any]] = {}
    pos = 0
    rest_chunks = []
    while True:
        m = re.compile(r"&(\w+)").search(joined, pos)
        if not m:
            rest_chunks.append(joined[pos:])
            break
        rest_chunks.append(joined[pos:m.start()])
        name = m.group(1).upper()
        # find the terminating '/' outside quotes
        q = None
        k = m.end()
        while k < len(joined):
            ch = joined[k]
            if q:
                if ch == q:
                    q = None
            elif ch in "'\"":
                q = ch
            elif ch == "/":
                break
            k += 1
        else:
            raise ConversionError("namelist &%s is not terminated by '/'" % name)
        namelists[name] = _parse_namelist_body(joined[m.end():k])
        pos = k + 1
    cards: Dict[str, Card] = {}
    current: Optional[Card] = None
    for line in "\n".join(rest_chunks).splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        mc = _CARD_RE.match(s)
        if mc:
            current = Card(mc.group(1).upper(), mc.group(2).lower())
            cards[current.name] = current
            continue
        if current is None:
            raise ConversionError("unexpected text outside a namelist or card: %r" % s)
        current.lines.append(s)
    return PwInput(namelists, cards)


_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.Pow: operator.pow, ast.USub: operator.neg,
        ast.UAdd: operator.pos}


def eval_number(text: str) -> float:
    """Evaluate a pw.x coordinate: a number or a simple expression such as ``1/3``."""
    t = re.sub(r"(\d|\.)[dD]([+-]?\d)", r"\1e\2", text.strip()).replace("^", "**")
    try:
        return float(t)
    except ValueError:
        pass

    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.BinOp) and type(node.op) in _OPS:
            return _OPS[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp) and type(node.op) in _OPS:
            return _OPS[type(node.op)](ev(node.operand))
        raise ValueError(text)

    try:
        return float(ev(ast.parse(t, mode="eval")))
    except (SyntaxError, ValueError, ZeroDivisionError):
        raise ConversionError("cannot read the number %r" % text)


# ---------------------------------------------------------------------------------------------
# Cell
# ---------------------------------------------------------------------------------------------


def _celldm(system: Dict[str, Any], ibrav: int) -> Dict[int, float]:
    """celldm(1..6) from either celldm(i) or A, B, C, cosAB, cosAC, cosBC (QE abc2celldm)."""
    cd = {i: float(system["celldm(%d)" % i]) for i in range(1, 7) if "celldm(%d)" % i in system}
    abc = {k: system[k] for k in ("a", "b", "c", "cosab", "cosac", "cosbc") if k in system}
    if cd and abc:
        raise ConversionError("celldm and A,B,C,cosAB,... are mutually exclusive")
    if abc:
        if "a" not in abc:
            raise ConversionError("A must be given when B, C or cosines are used")
        a = float(abc["a"]) * BOHR_PER_ANGSTROM
        cd = {1: a, 2: float(abc.get("b", 0.0)) * BOHR_PER_ANGSTROM / a,
              3: float(abc.get("c", 0.0)) * BOHR_PER_ANGSTROM / a}
        if ibrav == 14:
            cd[4], cd[5], cd[6] = (float(abc.get(k, 0.0)) for k in ("cosbc", "cosac", "cosab"))
        elif ibrav in (-12, -13):
            cd[5] = float(abc.get("cosac", 0.0))
        else:
            cd[4] = float(abc.get("cosab", 0.0))
    return cd


SUPPORTED_IBRAV = (1, 2, 3, -3, 4, 5, -5, 6, 7, 8, 9, -9, 10, 11, 12, -12, 14)


def qe_lattice(ibrav: int, cd: Dict[int, float]) -> np.ndarray:
    """Lattice vectors (rows, bohr) of QE's Bravais lattice ``ibrav`` (QE latgen.f90)."""
    if ibrav not in SUPPORTED_IBRAV:
        raise ConversionError("ibrav=%d is not supported by the converter; use ibrav=0 with "
                              "CELL_PARAMETERS" % ibrav)
    a = cd.get(1, 0.0)
    if a <= 0:
        raise ConversionError("celldm(1) (or A) must be positive for ibrav=%d" % ibrav)
    ba, ca = cd.get(2, 0.0), cd.get(3, 0.0)
    c4, c5, c6 = cd.get(4, 0.0), cd.get(5, 0.0), cd.get(6, 0.0)
    b, c = ba * a, ca * a

    def need(*vals):
        for name, v in vals:
            if v <= 0:
                raise ConversionError("%s must be positive for ibrav=%d" % (name, ibrav))

    if ibrav == 1:
        v = [[a, 0, 0], [0, a, 0], [0, 0, a]]
    elif ibrav == 2:
        h = a / 2
        v = [[-h, 0, h], [0, h, h], [-h, h, 0]]
    elif ibrav == 3:
        h = a / 2
        v = [[h, h, h], [-h, h, h], [-h, -h, h]]
    elif ibrav == -3:
        h = a / 2
        v = [[-h, h, h], [h, -h, h], [h, h, -h]]
    elif ibrav == 4:
        need(("celldm(3)", ca))
        v = [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, c]]
    elif ibrav in (5, -5):
        if not -0.5 < c4 < 1.0:
            raise ConversionError("celldm(4) = cos(gamma) must lie in (-0.5, 1) for ibrav=%d" % ibrav)
        tx = math.sqrt((1 - c4) / 2)
        ty = math.sqrt((1 - c4) / 6)
        tz = math.sqrt((1 + 2 * c4) / 3)
        if ibrav == 5:
            v = [[a * tx, -a * ty, a * tz], [0, 2 * a * ty, a * tz], [-a * tx, -a * ty, a * tz]]
        else:
            ap = a / math.sqrt(3)
            u = tz - 2 * math.sqrt(2) * ty
            w = tz + math.sqrt(2) * ty
            v = [[ap * u, ap * w, ap * w], [ap * w, ap * u, ap * w], [ap * w, ap * w, ap * u]]
    elif ibrav == 6:
        need(("celldm(3)", ca))
        v = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    elif ibrav == 7:
        need(("celldm(3)", ca))
        h = a / 2
        v = [[h, -h, c / 2], [h, h, c / 2], [-h, -h, c / 2]]
    elif ibrav == 8:
        need(("celldm(2)", ba), ("celldm(3)", ca))
        v = [[a, 0, 0], [0, b, 0], [0, 0, c]]
    elif ibrav == 9:
        need(("celldm(2)", ba), ("celldm(3)", ca))
        v = [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [0, 0, c]]
    elif ibrav == -9:
        need(("celldm(2)", ba), ("celldm(3)", ca))
        v = [[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c]]
    elif ibrav == 10:
        need(("celldm(2)", ba), ("celldm(3)", ca))
        v = [[a / 2, 0, c / 2], [a / 2, b / 2, 0], [0, b / 2, c / 2]]
    elif ibrav == 11:
        need(("celldm(2)", ba), ("celldm(3)", ca))
        v = [[a / 2, b / 2, c / 2], [-a / 2, b / 2, c / 2], [-a / 2, -b / 2, c / 2]]
    elif ibrav == 12:
        need(("celldm(2)", ba), ("celldm(3)", ca))
        s = math.sqrt(1 - c4 * c4)
        v = [[a, 0, 0], [b * c4, b * s, 0], [0, 0, c]]
    elif ibrav == -12:
        need(("celldm(2)", ba), ("celldm(3)", ca))
        s = math.sqrt(1 - c5 * c5)
        v = [[a, 0, 0], [0, b, 0], [c * c5, 0, c * s]]
    else:  # 14
        need(("celldm(2)", ba), ("celldm(3)", ca))
        cal, cbe, cga = c4, c5, c6
        sga = math.sqrt(1 - cga * cga)
        term = 1 + 2 * cal * cbe * cga - cal ** 2 - cbe ** 2 - cga ** 2
        if term <= 0 or sga == 0:
            raise ConversionError("celldm(4..6) do not describe a valid triclinic cell")
        v = [[a, 0, 0], [b * cga, b * sga, 0],
             [c * cbe, c * (cal - cbe * cga) / sga, c * math.sqrt(term) / sga]]
    return np.array(v, dtype=float)


def build_cell(pw: PwInput) -> Tuple[np.ndarray, float]:
    """Return (lattice in bohr, alat in bohr)."""
    system = pw.nl("SYSTEM")
    if "ibrav" not in system:
        raise ConversionError("ibrav is missing from &SYSTEM")
    ibrav = int(system["ibrav"])
    cd = _celldm(system, ibrav)
    card = pw.cards.get("CELL_PARAMETERS")
    if ibrav != 0:
        if card is not None:
            raise ConversionError("CELL_PARAMETERS given together with ibrav=%d" % ibrav)
        return qe_lattice(ibrav, cd), cd[1]
    if card is None:
        raise ConversionError("ibrav=0 needs a CELL_PARAMETERS card")
    if len(card.lines) < 3:
        raise ConversionError("CELL_PARAMETERS needs three lines")
    raw = np.array([[eval_number(x) for x in l.split()[:3]] for l in card.lines[:3]])
    opt = card.option or ("alat" if 1 in cd else "bohr")
    if opt == "alat":
        if 1 not in cd:
            raise ConversionError("CELL_PARAMETERS alat needs celldm(1) or A")
        lat = raw * cd[1]
    elif opt == "bohr":
        if cd:
            raise ConversionError("CELL_PARAMETERS bohr cannot be combined with celldm or A")
        lat = raw
    elif opt == "angstrom":
        if cd:
            raise ConversionError("CELL_PARAMETERS angstrom cannot be combined with celldm or A")
        lat = raw * BOHR_PER_ANGSTROM
    else:
        raise ConversionError("unknown CELL_PARAMETERS option %r" % opt)
    alat = cd.get(1, float(np.linalg.norm(lat[0])))
    return lat, alat


# ---------------------------------------------------------------------------------------------
# Keyword tables
# ---------------------------------------------------------------------------------------------

#: QE input_dft (lower case) -> RMG exchange_correlation_type (values of Input/InputOpts.cpp).
XC_MAP = {k: k for k in (
    "pz", "bp", "pw91", "blyp", "pbe", "revpbe", "pw86pbe", "b86bpbe", "pbesol", "q2d", "hcth",
    "olyp", "wc", "sogga", "optbk88", "optb86b", "ev93", "tpss", "m06l", "tb09", "pbe0", "hse",
    "b3lyp", "gaupbe", "vdw-df", "vdw-df-cx", "sla+pw+pbe+vdw1", "vdw-df-c09", "scan", "sca0",
    "r2scan")}
XC_MAP.update({"lda": "pz", "sla+pz": "pz", "sla+pw+pbx+pbc": "pbe", "hf": "hartree-fock",
               "hartree-fock": "hartree-fock", "pbe0": "pbe0", "hse06": "hse"})
HYBRIDS = {"pbe0", "hse", "b3lyp", "gaupbe", "hartree-fock", "sca0"}
LIBXC_ONLY = {"scan", "sca0", "r2scan"}
HSE_OMEGA_BOHR = 0.106

SMEARING_MAP = {
    "gaussian": ("MethfesselPaxton", 0), "gauss": ("MethfesselPaxton", 0),
    "methfessel-paxton": ("MethfesselPaxton", 1), "m-p": ("MethfesselPaxton", 1),
    "mp": ("MethfesselPaxton", 1),
    "marzari-vanderbilt": ("Cold Smearing", None), "cold": ("Cold Smearing", None),
    "m-v": ("Cold Smearing", None), "mv": ("Cold Smearing", None),
    "fermi-dirac": ("Fermi Dirac", None), "f-d": ("Fermi Dirac", None),
    "fd": ("Fermi Dirac", None),
}
TETRA_MAP = {"tetrahedra": "Bloechl", "tetrahedra_lin": "Linear", "tetrahedra_opt": "Optimized"}

CALC_MAP = {"scf": "Quench Electrons", "nscf": "NSCF", "relax": "Relax Structure",
            "vc-relax": "Relax Structure", "md": "Constant Volume And Energy"}

# Keys the converter understands, per namelist. Anything else is reported.
HANDLED = {
    "CONTROL": {"calculation", "title", "prefix", "outdir", "pseudo_dir", "verbosity",
                "restart_mode", "wf_collect", "disk_io", "tprnfor", "tstress", "etot_conv_thr",
                "forc_conv_thr", "nstep", "iprint", "max_seconds", "wfcdir", "lkpoint_dir"},
    "SYSTEM": {"ibrav", "nat", "ntyp", "ecutwfc", "ecutrho", "ecutfock", "nbnd", "tot_charge",
               "tot_magnetization", "nspin", "noncolin", "lspinorb", "occupations", "smearing",
               "degauss", "input_dft", "exx_fraction", "screening_parameter",
               "exxdiv_treatment", "x_gamma_extrapolation", "ecutvcut", "nqx1", "nqx2", "nqx3",
               "nosym", "noinv", "vdw_corr", "dftd3_version", "dftd3_threebody", "london",
               "a", "b", "c", "cosab", "cosac", "cosbc", "nr1", "nr2", "nr3", "nr1s", "nr2s",
               "nr3s", "ace", "starting_spin_angle", "use_all_frac", "force_symmorphic",
               "nosym_evc", "assume_isolated", "lda_plus_u", "constrained_magnetization",
               "lforcet"},
    "ELECTRONS": {"conv_thr", "electron_maxstep", "mixing_beta", "mixing_mode", "mixing_ndim",
                  "diagonalization", "startingwfc", "startingpot", "exx_maxstep",
                  "scf_must_converge", "diago_david_ndim", "diago_thr_init", "diago_full_acc",
                  "mixing_fixed_ns"},
    "IONS": {"ion_dynamics", "ion_temperature", "tempw", "upscale", "bfgs_ndim",
             "trust_radius_max", "trust_radius_min", "trust_radius_ini", "pot_extrapolation",
             "wfc_extrapolation"},
    "CELL": {"cell_dynamics", "cell_dofree", "press", "press_conv_thr", "wmass", "cell_factor"},
}
INDEXED_HANDLED = {"celldm", "starting_magnetization", "angle1", "angle2"}
RUN_CONTROL_ONLY = {"prefix", "outdir", "verbosity", "restart_mode", "wf_collect",
                    "disk_io", "iprint", "wfcdir", "lkpoint_dir", "nr1", "nr2", "nr3", "nr1s",
                    "nr2s", "nr3s", "ace", "diago_david_ndim", "diago_thr_init",
                    "diago_full_acc", "scf_must_converge", "upscale", "bfgs_ndim",
                    "trust_radius_max", "trust_radius_min", "trust_radius_ini",
                    "pot_extrapolation", "wfc_extrapolation", "use_all_frac",
                    "force_symmorphic", "nosym_evc", "lforcet", "startingwfc", "startingpot",
                    "diagonalization", "tprnfor"}
PHYSICS_UNSUPPORTED = {
    "CONTROL": {"tefield", "lelfield", "lberry", "gate", "lfcp", "trism", "dipfield", "lorbm",
                "nberrycyc", "gdir", "nppstr"},
    "SYSTEM": {"hubbard_u", "hubbard_j0", "hubbard_alpha", "hubbard_beta", "starting_charge",
               "one_atom_occupations", "fixed_magnetization", "b_field", "lambda", "report",
               "edir", "emaxpos", "eopreg", "eamp", "qcutz", "q2sigma", "ecfixed", "xdm",
               "xdm_a1", "xdm_a2", "ts_vdw_econv_thr", "ts_vdw_isolated", "esm_bc",
               "esm_w", "esm_efield", "esm_nfit", "lgcscf", "block", "zgate", "relaxz",
               "ensemble_energies", "space_group", "uniqueb", "origin_choice",
               "rhombohedral", "nextffield", "sic_gamma", "pol_type", "sci_vb", "sci_cb",
               "lspinorb_sic", "london_s6", "london_rcut", "london_c6", "london_rvdw",
               "gcscf_mu", "gcscf_conv_thr", "gcscf_beta", "hubbard_occ"},
    "ELECTRONS": {"efield", "efield_cart", "efield_phase", "tqr", "real_space",
                  "occupation_constraints"},
}


# ---------------------------------------------------------------------------------------------
# Converter
# ---------------------------------------------------------------------------------------------


def _base_key(k: str) -> str:
    return re.sub(r"\(.*\)$", "", k)


def _review_keys(pw: PwInput, notes: Notes) -> None:
    for name, nl in pw.namelists.items():
        if name not in ("CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL"):
            notes.unsupported("namelist &%s is not converted" % name)
            continue
        for key in nl:
            base = _base_key(key)
            if base in PHYSICS_UNSUPPORTED.get(name, set()):
                notes.unsupported("&%s %s is not converted (no RMG equivalent in this converter)"
                                  % (name, key))
            elif key in HANDLED.get(name, set()) or base in INDEXED_HANDLED:
                if base in RUN_CONTROL_ONLY:
                    notes.info("&%s %s is run control only; not converted" % (name, key))
            else:
                notes.warn("&%s %s is not recognised by the converter and was ignored"
                           % (name, key))
    for cname in ("OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES", "SOLVENTS", "HUBBARD",
                  "ADDITIONAL_K_POINTS", "TOTAL_CHARGE"):
        if cname in pw.cards:
            notes.unsupported("card %s is not converted" % cname)
    if "ATOMIC_VELOCITIES" in pw.cards:
        notes.warn("ATOMIC_VELOCITIES not converted; RMG draws initial velocities itself")


def _read_species(pw: PwInput, pseudo_search: Optional[str], notes: Notes):
    card = pw.cards.get("ATOMIC_SPECIES")
    if card is None:
        raise ConversionError("ATOMIC_SPECIES card is missing")
    ntyp = int(pw.nl("SYSTEM").get("ntyp", len(card.lines)))
    if len(card.lines) < ntyp:
        raise ConversionError("ATOMIC_SPECIES has %d lines, ntyp = %d" % (len(card.lines), ntyp))
    species = []
    for line in card.lines[:ntyp]:
        f = line.split()
        if len(f) < 3:
            raise ConversionError("bad ATOMIC_SPECIES line %r" % line)
        label, upf = f[0], f[2]
        header: Optional[UpfHeader] = None
        if pseudo_search is not None:
            path = os.path.join(pseudo_search, upf)
            if os.path.isfile(path):
                try:
                    header = read_upf_header(path)
                except (OSError, ValueError) as exc:
                    notes.warn("could not read %s: %s" % (path, exc))
            else:
                notes.warn("pseudopotential %s not found; its valence charge and format could "
                           "not be checked" % path)
        elem = (header.element if header and header.element else None) or element_from_label(label)
        if elem is None or element_from_label(elem) != elem:
            raise ConversionError("cannot tell the element of species %r" % label)
        if header is not None:
            check_upf_for_rmg(header, notes)
        species.append({"label": label, "element": elem, "upf": upf, "header": header})
    return species


def _positions(pw: PwInput, lattice: np.ndarray, alat: float, labels: List[str]):
    card = pw.cards.get("ATOMIC_POSITIONS")
    if card is None:
        raise ConversionError("ATOMIC_POSITIONS card is missing")
    nat = int(pw.nl("SYSTEM").get("nat", len(card.lines)))
    if len(card.lines) < nat:
        raise ConversionError("ATOMIC_POSITIONS has %d lines, nat = %d" % (len(card.lines), nat))
    opt = card.option or "alat"
    names, xyz, flags = [], [], []
    for line in card.lines[:nat]:
        f = line.split()
        if f[0] not in labels:
            raise ConversionError("atom label %r is not in ATOMIC_SPECIES" % f[0])
        names.append(f[0])
        xyz.append([eval_number(x) for x in f[1:4]])
        flags.append([int(x) for x in f[4:7]] if len(f) >= 7 else [1, 1, 1])
    xyz = np.array(xyz, dtype=float)
    inv = np.linalg.inv(lattice)
    if opt == "crystal":
        frac = xyz
    elif opt == "alat":
        frac = (xyz * alat) @ inv
    elif opt == "bohr":
        frac = xyz @ inv
    elif opt == "angstrom":
        frac = (xyz * BOHR_PER_ANGSTROM) @ inv
    else:
        raise ConversionError("ATOMIC_POSITIONS {%s} is not supported" % opt)
    return names, frac, np.array(flags, dtype=int)


def _kpoints(pw: PwInput, lattice: np.ndarray, alat: float, notes: Notes):
    """Return ('mesh', mesh, shift) or ('list', kfrac, weights)."""
    card = pw.cards.get("K_POINTS")
    if card is None:
        raise ConversionError("K_POINTS card is missing")
    opt = card.option or "tpiba"
    if opt == "gamma":
        return ("mesh", (1, 1, 1), (0, 0, 0))
    if opt == "automatic":
        if not card.lines:
            raise ConversionError("K_POINTS automatic needs a line nk1 nk2 nk3 sk1 sk2 sk3")
        f = [int(x) for x in card.lines[0].split()[:6]]
        if len(f) != 6 or min(f[:3]) < 1 or any(s not in (0, 1) for s in f[3:]):
            raise ConversionError("bad K_POINTS automatic line %r" % card.lines[0])
        return ("mesh", tuple(f[:3]), tuple(f[3:]))
    if opt in ("tpiba", "crystal"):
        n = int(card.lines[0].split()[0])
        rows = [[eval_number(x) for x in l.split()[:4]] for l in card.lines[1:1 + n]]
        if len(rows) != n or any(len(r) != 4 for r in rows):
            raise ConversionError("K_POINTS %s list is incomplete" % opt)
        k = np.array([r[:3] for r in rows])
        w = [r[3] for r in rows]
        if opt == "tpiba":
            k = k @ lattice.T / alat  # k_frac_i = k . a_i / alat for k in 2 pi/alat
        notes.info("explicit k-point list written in reciprocal-lattice units; RMG normalises "
                   "the weights and does not symmetrise the list")
        return ("list", k, w)
    raise ConversionError("K_POINTS %s is not supported (band paths need a manual RMG "
                          "kpoints_bandstructure block)" % opt)


def convert(text: str, input_dir: str = ".", pseudo_dir: Optional[str] = None,
            grid: Optional[Sequence[int]] = None, grid_spacing: Optional[float] = None,
            grid_multiple: int = 4, lenient: bool = False) -> Tuple[RmgDeck, Notes]:
    """Convert pw.x input text. Returns the RMG deck and the conversion notes.

    ``input_dir`` resolves a relative ``pseudo_dir``; ``pseudo_dir`` overrides the one in the
    input (and is written to the deck).
    """
    notes = Notes(lenient)
    pw = parse_pw_input(text)
    control, system, electrons = pw.nl("CONTROL"), pw.nl("SYSTEM"), pw.nl("ELECTRONS")
    ions, cellnl = pw.nl("IONS"), pw.nl("CELL")
    _review_keys(pw, notes)
    deck = RmgDeck()

    # -- description ---------------------------------------------------------------------
    title = control.get("title") or control.get("prefix") or ""
    deck.set("description", sanitize_description("converted from QE pw.x input " + str(title)))

    # -- calculation ---------------------------------------------------------------------
    calc = str(control.get("calculation", "scf")).lower()
    if calc not in CALC_MAP:
        notes.unsupported("calculation = '%s' has no RMG counterpart in this converter" % calc)
        calc = "scf"
    mode = CALC_MAP[calc]
    if calc == "md":
        tcontrol = str(ions.get("ion_temperature", "not_controlled")).lower()
        if tcontrol == "nose":
            mode = "Constant Temperature And Energy"
            deck.set("md_temperature_control", "Nose Hoover Chains")
            if "tempw" in ions:
                deck.set("md_temperature", float(ions["tempw"]))
        elif tcontrol != "not_controlled":
            notes.unsupported("ion_temperature = '%s' is not converted" % tcontrol)
        notes.warn("the MD time step dt is not converted; set ionic_time_step for RMG by hand")
    deck.set("calculation_mode", mode)
    if calc == "vc-relax":
        dof = str(cellnl.get("cell_dofree", "all")).lower()
        if dof != "all":
            notes.unsupported("cell_dofree = '%s' is not converted (RMG cell_movable needs a "
                              "manual 3x3 mask)" % dof)
        deck.set("cell_relax", True)
        if "press" in cellnl and float(cellnl["press"]) != 0.0:
            notes.unsupported("a target pressure (press) is not converted")
    if calc in ("relax", "vc-relax", "md"):
        deck.set("max_md_steps", int(control.get("nstep", 50)))
    if calc in ("relax", "vc-relax"):
        deck.set("relax_max_force", float(control.get("forc_conv_thr", 1.0e-3)) / 2.0)
        dyn = str(ions.get("ion_dynamics", "bfgs")).lower()
        rm = {"bfgs": "LBFGS", "fire": "FIRE"}.get(dyn)
        if rm:
            deck.set("relax_method", rm)
        else:
            notes.warn("ion_dynamics = '%s' has no direct RMG match; RMG's default relax_method "
                       "(LBFGS) is used" % dyn)
        if "etot_conv_thr" in control:
            notes.info("etot_conv_thr not converted; RMG stops relaxations on relax_max_force")
    if control.get("tstress"):
        deck.set("stress", True)
    if "max_seconds" in control:
        deck.set("max_run_seconds", float(control["max_seconds"]))
        notes.info("max_run_seconds is a recent RMG key; older RMG builds ignore it")

    # -- cell, species, atoms --------------------------------------------------------------
    lattice, alat = build_cell(pw)
    ps_dir_in = pseudo_dir if pseudo_dir is not None else control.get("pseudo_dir")
    if ps_dir_in is None:
        ps_dir_in = os.environ.get("ESPRESSO_PSEUDO") or os.path.join(
            os.path.expanduser("~"), "espresso", "pseudo")
        notes.warn("pseudo_dir not set; using QE's default %s" % ps_dir_in)
    ps_search = ps_dir_in if os.path.isabs(ps_dir_in) else os.path.join(input_dir, ps_dir_in)
    species = _read_species(pw, ps_search, notes)
    labels = [s["label"] for s in species]
    by_label = {s["label"]: s for s in species}
    names, frac, flags = _positions(pw, lattice, alat, labels)
    if calc in ("scf", "nscf") and (flags == 0).any():
        notes.info("if_pos constraints carried into RMG movable flags (inactive for %s)" % calc)
    elif (flags == 0).any():
        notes.warn("if_pos constraints carried into RMG per-axis movable flags; check them")

    # pseudopotentials: one file per element
    pp: Dict[str, str] = {}
    for s in species:
        if s["element"] in pp and pp[s["element"]] != s["upf"]:
            notes.unsupported("species of element %s use different pseudopotentials (%s, %s); "
                              "RMG keys pseudopotentials by element"
                              % (s["element"], pp[s["element"]], s["upf"]))
        pp.setdefault(s["element"], s["upf"])
    if len({s["element"] for s in species}) != len(species):
        notes.info("several QE species map to one element; their differences (e.g. starting "
                   "magnetization) are kept per atom")

    # -- spin ------------------------------------------------------------------------------
    nspin = int(system.get("nspin", 1))
    noncolin = bool(system.get("noncolin", False)) or bool(system.get("lspinorb", False))
    if bool(system.get("lspinorb", False)) and not bool(system.get("noncolin", False)):
        notes.warn("lspinorb without noncolin: treated as a noncollinear spin-orbit run")
    if nspin not in (1, 2):
        raise ConversionError("nspin = %d is not valid" % nspin)
    if noncolin and nspin == 2:
        raise ConversionError("noncolin and nspin=2 are mutually exclusive")
    if str(system.get("constrained_magnetization", "none")).lower() != "none":
        notes.unsupported("constrained_magnetization is not converted")
    if system.get("lda_plus_u"):
        notes.unsupported("DFT+U (lda_plus_u / HUBBARD) is not converted; RMG's Hubbard_U "
                          "input needs orbital labels and differs in its projectors")
    if str(system.get("assume_isolated", "none")).lower() not in ("none",):
        notes.unsupported("assume_isolated = '%s' is not converted" % system["assume_isolated"])

    spin_cols = 0
    spins = [()] * len(names)
    if nspin == 2 or noncolin:
        spin_cols = 3 if noncolin else 1
        spins = []
        for nm in names:
            it = labels.index(nm) + 1
            m = float(system.get("starting_magnetization(%d)" % it, 0.0))
            r = clamp_spin_ratio(m / 2.0, notes, "species " + nm)
            if noncolin:
                spins.append((r, float(system.get("angle1(%d)" % it, 0.0)),
                              float(system.get("angle2(%d)" % it, 0.0))))
            else:
                spins.append((r,))
        if noncolin:
            deck.set("noncollinear", True)
            if bool(system.get("lspinorb", False)):
                deck.set("spinorbit", True)
                for s in species:
                    h = s["header"]
                    if h is not None and h.has_so is False:
                        notes.warn("%s is not fully relativistic (has_so false); RMG needs "
                                   "fully relativistic pseudopotentials for spin-orbit"
                                   % s["upf"])
        else:
            deck.set("spin_polarization", True)
        notes.info("RMG initial spin column = starting_magnetization / 2, i.e. the same "
                   "up/down split of the atomic density")

    elements = [by_label[n]["element"] for n in names]
    st = Structure(lattice, elements, frac, flags, spins)
    st, orient = apply_orientation(st, notes)

    # -- electrons count -------------------------------------------------------------------
    tot_charge = float(system.get("tot_charge", 0.0))
    zvals = [by_label[n]["header"].z_valence if by_label[n]["header"] else None for n in names]
    nel = None if any(z is None for z in zvals) else sum(zvals) - tot_charge
    if tot_charge != 0.0:
        deck.set("system_charge", tot_charge)

    # -- functional ------------------------------------------------------------------------
    xc = None
    if "input_dft" in system:
        key = str(system["input_dft"]).strip().lower()
        if key not in XC_MAP:
            notes.unsupported("input_dft = '%s' has no RMG exchange_correlation_type" % key)
        else:
            xc = XC_MAP[key]
            deck.set("exchange_correlation_type", xc)
            if xc in LIBXC_ONLY:
                notes.warn("exchange_correlation_type %s needs an RMG build with libxc" % xc)
    else:
        notes.info("input_dft not set: RMG, like QE, takes the functional from the "
                   "pseudopotentials (AUTO_XC)")
    hybrid = xc in HYBRIDS
    if "exx_fraction" in system:
        deck.set("exx_fraction", float(system["exx_fraction"]))
    if "screening_parameter" in system:
        om = float(system["screening_parameter"])
        if xc == "hse" and abs(om - HSE_OMEGA_BOHR) > 1e-9:
            notes.unsupported("screening_parameter = %g bohr^-1: RMG has no input for the HSE "
                              "screening length and always uses %g" % (om, HSE_OMEGA_BOHR))
        elif xc != "hse":
            notes.warn("screening_parameter ignored for input_dft '%s'" % xc)
    if "ecutfock" in system:
        notes.info("ecutfock has no RMG equivalent; not converted")

    # -- occupations -----------------------------------------------------------------------
    occ = str(system.get("occupations", "fixed")).lower()
    smeared = False
    if occ == "smearing":
        sm = str(system.get("smearing", "gaussian")).lower()
        if sm not in SMEARING_MAP:
            raise ConversionError("unknown smearing '%s'" % sm)
        otype, order = SMEARING_MAP[sm]
        deck.set("occupations_type", otype)
        if order is not None:
            deck.set("MP_order", order)
        width_ev = float(system.get("degauss", 0.0)) * RY_EV
        if width_ev > 2.0:
            notes.unsupported("degauss = %.4g eV exceeds RMG's 2 eV limit" % width_ev)
        deck.set("occupation_electron_temperature_eV", width_ev)
        if sm in ("gaussian", "gauss"):
            notes.info("QE Gaussian smearing = RMG MethfesselPaxton with MP_order 0 "
                       "(RMG's 'Gaussian' occupation is a different function)")
        smeared = True
    elif occ in TETRA_MAP:
        deck.set("occupations_type", "Tetrahedron")
        deck.set("tetra_method", TETRA_MAP[occ])
        smeared = True
    elif occ == "fixed":
        deck.set("occupations_type", "Fixed")
    else:
        notes.unsupported("occupations = '%s' is not converted" % occ)

    tot_mag = system.get("tot_magnetization")
    if tot_mag is not None and float(tot_mag) == -10000:
        tot_mag = None
    nbnd = system.get("nbnd")
    unocc = None
    if nbnd is not None:
        if nel is None:
            notes.warn("nbnd = %s not converted: the valence charge of every pseudopotential "
                       "is needed (UPF files not found); RMG's default number of unoccupied "
                       "states is used" % nbnd)
        else:
            nocc = nocc_rmg(nel, noncolin)
            unocc = int(nbnd) - nocc
            if unocc < 0:
                notes.warn("nbnd = %d is below RMG's %d occupied states; unoccupied states set "
                           "to 0" % (int(nbnd), nocc))
                unocc = 0
            deck.set("unoccupied_states_per_kpoint", unocc)
    if tot_mag is not None:
        if nspin != 2:
            raise ConversionError("tot_magnetization requires nspin = 2")
        tot_mag = float(tot_mag)
        if smeared:
            deck.set("tot_magnetization", tot_mag)
            notes.info("tot_magnetization is a recent RMG key; older RMG builds ignore it silently")
            if nel is not None and unocc is not None:
                if (nel + abs(tot_mag)) / 2.0 > nocc_rmg(nel) + unocc:
                    notes.warn("nbnd leaves too few states for the majority spin at "
                               "tot_magnetization = %g" % tot_mag)
        else:
            if nel is None:
                notes.unsupported("fixed occupations with tot_magnetization need explicit RMG "
                                  "spin occupations, which need the valence charges (UPF files "
                                  "not found)")
            else:
                nup, ndn = (nel + tot_mag) / 2.0, (nel - tot_mag) / 2.0
                if abs(nup - round(nup)) > 1e-8 or abs(ndn - round(ndn)) > 1e-8:
                    raise ConversionError("fixed occupations need integer spin populations; "
                                          "got %g up and %g down" % (nup, ndn))
                nup, ndn = int(round(nup)), int(round(ndn))
                nst = max(nup, ndn) if unocc is None else max(nup, ndn, nocc_rmg(nel) + unocc)
                deck.set("states_count_and_occupation_spin_up",
                         "%d 1.0 %d 0.0" % (nup, nst - nup) if nst > nup else "%d 1.0" % nup)
                deck.set("states_count_and_occupation_spin_down",
                         "%d 1.0 %d 0.0" % (ndn, nst - ndn) if nst > ndn else "%d 1.0" % ndn)
    elif nspin == 2 and occ == "fixed":
        raise ConversionError("nspin = 2 with fixed occupations needs tot_magnetization (as in "
                              "pw.x)")

    # -- k-points --------------------------------------------------------------------------
    kp = _kpoints(pw, lattice, alat, notes)
    if kp[0] == "mesh":
        kmesh, kshift = kp[1], kp[2]
        deck.set("kpoint_mesh", list(kmesh))
        deck.set("kpoint_is_shift", list(kshift))
    else:
        kmesh = None
        kfrac = np.asarray(kp[1]) @ orient.basis_change.T
        write_kpoint_list(deck, kfrac, kp[2])

    # -- exact exchange --------------------------------------------------------------------
    nq = [int(system.get("nqx%d" % i, 1)) for i in (1, 2, 3)]
    if hybrid:
        if kmesh is None:
            notes.unsupported("hybrid functionals with an explicit k-point list are not "
                              "converted (RMG exact exchange needs a k mesh)")
        else:
            check_exx_qmesh(kmesh, nq)
            deck.set("exx_qpoint_mesh", nq)
            notes.info("exx_qpoint_mesh is written explicitly (QE default nqx = 1 1 1; RMG's "
                       "default is the full k mesh); it needs an RMG build with "
                       "exx_qpoint_mesh support")
        deck.set("x_gamma_extrapolation", bool(system.get("x_gamma_extrapolation", True)))
        div = str(system.get("exxdiv_treatment", "gygi-baldereschi")).lower()
        if div not in ("gygi-baldereschi", "none"):
            notes.unsupported("exxdiv_treatment = '%s' is not available in RMG" % div)
        else:
            deck.set("exxdiv_treatment", div)
        if "exx_maxstep" in electrons:
            deck.set("max_exx_steps", int(electrons["exx_maxstep"]))
        if "ecutvcut" in system:
            notes.info("ecutvcut not converted")
    elif any(n != 1 for n in nq):
        notes.warn("nqx1..3 ignored: the functional is not a hybrid")

    # -- vdW -------------------------------------------------------------------------------
    vdw = str(system.get("vdw_corr", "none")).lower()
    if system.get("london"):
        vdw = "grimme-d2"
    if vdw in ("grimme-d2", "dft-d", "d2"):
        deck.set("vdw_corr", "DFT-D2")
    elif vdw in ("grimme-d3", "dft-d3", "d3"):
        deck.set("vdw_corr", "DFT-D3")
        deck.set("dftd3_version", int(system.get("dftd3_version", 3)))
        if system.get("dftd3_threebody", True):
            notes.info("dftd3_threebody has no RMG input; RMG's DFT-D3 default is used")
    elif vdw != "none":
        notes.unsupported("vdw_corr = '%s' is not converted" % vdw)

    # -- symmetry --------------------------------------------------------------------------
    if system.get("nosym"):
        deck.set("use_symmetry", 0)
    if system.get("noinv"):
        deck.set("time_reversal", False)

    # -- grid ------------------------------------------------------------------------------
    if "ecutwfc" not in system and grid is None and grid_spacing is None:
        raise ConversionError("ecutwfc is missing")
    ecut = float(system["ecutwfc"]) if "ecutwfc" in system else None
    g = wavefunction_grid(st.lattice, ecut, grid_spacing, grid, grid_multiple, notes)
    deck.set("wavefunction_grid", list(g.grid))
    notes.info("wavefunction_grid %s %s; spacings %.4f %.4f %.4f bohr; equivalent cutoff "
               "%.2f Ry" % (" ".join(map(str, g.grid)), g.how, *g.spacing_bohr, g.ecut_equiv_ry))
    if ecut is not None:
        ecutrho = float(system.get("ecutrho", 4.0 * ecut))
        r = potential_grid_refinement(ecutrho, ecut)
        deck.set("potential_grid_refinement", r)
        if abs(r * r * ecut - ecutrho) > 1e-6 * ecutrho:
            notes.info("ecutrho = %g Ry -> potential_grid_refinement %d (density cutoff "
                       "%g x ecutwfc)" % (ecutrho, r, r * r))

    # -- SCF control -----------------------------------------------------------------------
    set_scf_thresholds(deck, float(electrons.get("conv_thr", 1.0e-6)) / 2.0, notes,
                       hybrid=hybrid, source="conv_thr", spin=(nspin == 2 or noncolin))
    deck.set("max_scf_steps", int(electrons.get("electron_maxstep", 100)))
    mix_mode = str(electrons.get("mixing_mode", "plain")).lower()
    deck.set("charge_mixing_type", "Broyden")
    if mix_mode != "plain":
        notes.warn("mixing_mode = '%s' mapped to RMG Broyden mixing" % mix_mode)
    deck.set("charge_density_mixing", float(electrons.get("mixing_beta", 0.7)))
    if "mixing_ndim" in electrons:
        deck.set("charge_broyden_order", max(1, min(12, int(electrons["mixing_ndim"]))))

    # -- pseudopotentials and structure ----------------------------------------------------
    deck.set("pseudo_dir", ps_dir_in)
    deck.set_block("pseudopotential", ["  %s  %s" % (e, f) for e, f in pp.items()])
    write_structure(deck, st, spin_cols)
    notes.info("atomic masses from ATOMIC_SPECIES are not converted (RMG uses its own table)")
    return deck, notes


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog=PROG, description=__doc__.split("\n\n")[0])
    ap.add_argument("input", help="pw.x input file")
    ap.add_argument("-o", "--output", default="-", help="RMG input file to write (default stdout)")
    ap.add_argument("--notes", default=None,
                    help="notes file (default <output>.notes; none when writing to stdout)")
    ap.add_argument("--pseudo-dir", default=None,
                    help="directory with the UPF files; overrides pseudo_dir of the input and is "
                         "written to the deck")
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument("--grid", nargs=3, type=int, metavar=("NX", "NY", "NZ"),
                     help="wavefunction grid to use instead of the one derived from ecutwfc")
    grp.add_argument("--grid-spacing", type=float, metavar="H",
                     help="target grid spacing in bohr instead of the ecutwfc mapping")
    ap.add_argument("--grid-multiple", type=int, default=4,
                    help="every grid dimension is a multiple of this (default 4)")
    ap.add_argument("--lenient", action="store_true",
                    help="drop unsupported features with a warning instead of stopping")
    a = ap.parse_args(argv)
    try:
        with open(a.input) as f:
            text = f.read()
        deck, notes = convert(text, os.path.dirname(os.path.abspath(a.input)), a.pseudo_dir,
                              a.grid, a.grid_spacing, a.grid_multiple, a.lenient)
        write_outputs(PROG, deck, notes, a.output, a.notes,
                      "%s notes for %s (RMG decks cannot hold comments)" % (PROG, a.input))
    except ConversionError as exc:
        print("%s: error: %s" % (PROG, exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
