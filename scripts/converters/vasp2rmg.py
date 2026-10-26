#!/usr/bin/env python3
"""Convert VASP inputs (POSCAR/CONTCAR, INCAR, KPOINTS) into an RMG (rmg-dft) input file.

Usage::

    vasp2rmg.py [-d DIR] [--poscar FILE] [--incar FILE] [--kpoints FILE] [--potcar FILE]
                [--species "Fe Si"] [--pseudo El=file.upf ...] [--system-charge Q]
                [--grid NX NY NZ | --grid-spacing H] [-o rmg.in] [--notes FILE] [--lenient]

By default the files POSCAR (else CONTCAR), INCAR, KPOINTS and POTCAR are taken from DIR
(default: the current directory); INCAR, KPOINTS and POTCAR are optional.

POTCAR files hold PAW datasets, which RMG cannot use. The POTCAR is only read for what does not
depend on the dataset: the species order (TITEL), the valence charges ZVAL (to turn NELECT into
a net charge) and ENMAX (VASP's default ENCUT). Without ``--pseudo`` the deck has no
pseudopotential block, so RMG uses its internal norm-conserving pseudopotentials; with
``--pseudo El=file.upf`` (one per element) the deck points RMG at those UPF files.

The deck goes to ``-o`` (or stdout). Everything the converter changed, dropped or could not
check is printed on stderr and written to ``<output>.notes``; RMG decks cannot carry comment
lines. Features that would change the physics and that RMG cannot represent make the
converter stop with exit status 2, unless ``--lenient`` is given.

Library use::

    from vasp2rmg import convert
    deck, notes = convert(open("POSCAR").read(), open("INCAR").read(), open("KPOINTS").read())
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rmg_input import (BOHR_PER_ANGSTROM, HA_EV, RY_EV, ConversionError, Notes,  # noqa: E402
                       RmgDeck, Structure, apply_orientation, check_exx_qmesh,
                       check_upf_for_rmg, clamp_spin_ratio, element_from_label, nocc_rmg,
                       read_upf_header, reciprocal, sanitize_description, set_scf_thresholds,
                       wavefunction_grid, write_kpoint_list, write_outputs, write_structure)

PROG = "vasp2rmg"

# ---------------------------------------------------------------------------------------------
# POSCAR
# ---------------------------------------------------------------------------------------------


@dataclass
class Poscar:
    comment: str
    lattice: np.ndarray            # Angstrom, rows
    symbols: List[str]             # one chemical symbol per atom
    frac: np.ndarray
    selective: Optional[np.ndarray]  # (N, 3) 0/1 or None


def _floats(line: str, n: int, what: str) -> List[float]:
    f = line.split()
    try:
        if len(f) >= n:
            return [float(x) for x in f[:n]]
    except ValueError:
        pass
    raise ConversionError("cannot read %s from POSCAR line %r" % (what, line))


def parse_poscar(text: str, species: Optional[Sequence[str]] = None,
                 fallback_species: Optional[Sequence[str]] = None) -> Poscar:
    """Parse a VASP 5 (species line) or VASP 4 POSCAR/CONTCAR.

    ``species`` (user given) must agree with a species line if there is one;
    ``fallback_species`` (e.g. from the POTCAR) is used only for a VASP 4 file.
    """
    lines = [l.rstrip() for l in text.splitlines()]
    while lines and not lines[-1].strip():
        lines.pop()
    if len(lines) < 8:
        raise ConversionError("POSCAR is too short")
    comment = lines[0].strip()
    scale = [float(x) for x in lines[1].split()[:3] if re.match(r"^[-+0-9.eEdD]+$", x)]
    if len(scale) not in (1, 3):
        raise ConversionError("bad POSCAR scaling line %r" % lines[1])
    lat = np.array([_floats(lines[i], 3, "a lattice vector") for i in (2, 3, 4)])
    i = 5
    tok = lines[i].split()
    if tok and not re.match(r"^\d+$", tok[0]):
        names = tok
        i += 1
        if species:
            if [element_from_label(s) for s in species] != [element_from_label(s) for s in names]:
                raise ConversionError("--species %s disagrees with the POSCAR species line %s"
                                      % (list(species), names))
    else:
        if not species and not fallback_species:
            raise ConversionError("VASP 4 POSCAR without a species line: give --species (or a "
                                  "POTCAR)")
        names = list(species or fallback_species)
    try:
        counts = [int(x) for x in lines[i].split()]
    except ValueError:
        raise ConversionError("bad POSCAR atom-count line %r" % lines[i])
    if len(counts) != len(names):
        raise ConversionError("POSCAR has %d species names and %d counts"
                              % (len(names), len(counts)))
    i += 1
    selective = False
    if lines[i].strip()[:1] in ("s", "S"):
        selective = True
        i += 1
    mode = lines[i].strip()[:1].lower()
    cartesian = mode in ("c", "k")
    i += 1
    n = sum(counts)
    if len(lines) < i + n:
        raise ConversionError("POSCAR lists %d atoms but has fewer coordinate lines" % n)
    xyz, flags = [], []
    for l in lines[i:i + n]:
        f = l.split()
        xyz.append(_floats(l, 3, "atomic coordinates"))
        if selective:
            if len(f) < 6 or any(x[:1].upper() not in "TF" for x in f[3:6]):
                raise ConversionError("selective dynamics flags missing on %r" % l)
            flags.append([1 if x[:1].upper() == "T" else 0 for x in f[3:6]])
    xyz = np.array(xyz)
    # scaling
    if len(scale) == 3:
        if min(scale) <= 0:
            raise ConversionError("three POSCAR scaling factors must all be positive")
        smat = np.diag(scale)
        lat = lat @ smat
        if cartesian:
            xyz = xyz @ smat
    else:
        s = scale[0]
        if s < 0:
            s = (abs(s) / abs(np.linalg.det(lat))) ** (1.0 / 3.0)
        elif s == 0:
            raise ConversionError("POSCAR scaling factor is zero")
        lat = lat * s
        if cartesian:
            xyz = xyz * s
    frac = xyz @ np.linalg.inv(lat) if cartesian else xyz
    symbols = []
    for nm, c in zip(names, counts):
        el = element_from_label(nm)
        if el is None:
            raise ConversionError("cannot tell the element of POSCAR species %r" % nm)
        symbols += [el] * c
    return Poscar(comment, lat, symbols, frac, np.array(flags, dtype=int) if selective else None)


# ---------------------------------------------------------------------------------------------
# INCAR
# ---------------------------------------------------------------------------------------------


def parse_incar(text: str) -> Dict[str, str]:
    """``{TAG: raw value}`` from INCAR text (comments ``#``/``!``, ``;`` and ``\\`` handled)."""
    out: Dict[str, str] = {}
    joined = re.sub(r"\\\s*\n", " ", text)
    for raw in joined.splitlines():
        line = re.split(r"[#!]", raw, 1)[0]
        for stmt in line.split(";"):
            if "=" not in stmt:
                if stmt.strip():
                    raise ConversionError("cannot read INCAR statement %r" % stmt.strip())
                continue
            k, v = stmt.split("=", 1)
            k = k.strip().upper()
            v = v.strip().strip('"').strip("'").strip()
            if not k:
                raise ConversionError("INCAR statement without a tag: %r" % stmt)
            out[k] = v
    return out


def incar_bool(v: str) -> bool:
    t = v.strip().strip(".").upper()
    if t in ("T", "TRUE"):
        return True
    if t in ("F", "FALSE"):
        return False
    raise ConversionError("not a logical value: %r" % v)


def incar_float(v: str) -> float:
    try:
        return float(v.split()[0].replace("d", "e").replace("D", "E"))
    except (ValueError, IndexError):
        raise ConversionError("not a number: %r" % v)


def incar_int(v: str) -> int:
    x = incar_float(v)
    if x != int(x):
        raise ConversionError("not an integer: %r" % v)
    return int(x)


def expand_magmom(v: str) -> List[float]:
    """Expand VASP's ``N*value`` shorthand: ``"2*3.0 1*-1 0"`` -> [3.0, 3.0, -1.0, 0.0]."""
    out: List[float] = []
    for tok in v.split():
        if "*" in tok:
            n, val = tok.split("*", 1)
            try:
                out += [float(val)] * int(n)
            except ValueError:
                raise ConversionError("bad MAGMOM entry %r" % tok)
        else:
            out.append(incar_float(tok))
    return out


# ---------------------------------------------------------------------------------------------
# KPOINTS and POTCAR
# ---------------------------------------------------------------------------------------------


def _shift_flags(mesh: Sequence[int], shift: Sequence[float], monkhorst: bool) -> Tuple[int, ...]:
    """RMG kpoint_is_shift for a VASP grid.

    RMG puts k-points at (i + 0.5 s)/N with s in {0, 1}. A Gamma-centred VASP grid with an
    extra shift t (fraction of a grid step) needs t in {0, 1/2}. A Monkhorst-Pack grid is
    Gamma-centred for odd N and shifted by half a step for even N, plus any extra shift.
    """
    flags = []
    for n, t in zip(mesh, shift):
        total = (t + (0.5 if (monkhorst and n % 2 == 0) else 0.0)) % 1.0
        if abs(total) < 1e-8 or abs(total - 1.0) < 1e-8:
            flags.append(0)
        elif abs(total - 0.5) < 1e-8:
            flags.append(1)
        else:
            raise ConversionError("k-point shift %g of a grid step is not representable in RMG "
                                  "(only 0 or 1/2)" % total)
    return tuple(flags)


def mesh_from_length(lattice_ang: np.ndarray, rk: float) -> Tuple[int, int, int]:
    """VASP fully automatic mode: N_i = max(1, int(R_k |b_i| + 0.5)), b without 2 pi."""
    b = np.linalg.norm(reciprocal(lattice_ang), axis=1)
    return tuple(max(1, int(rk * x + 0.5)) for x in b)


def mesh_from_kspacing(lattice_ang: np.ndarray, kspacing: float) -> Tuple[int, int, int]:
    """INCAR KSPACING: N_i = max(1, ceil(2 pi |b_i| / KSPACING)), b without 2 pi."""
    b = np.linalg.norm(reciprocal(lattice_ang), axis=1)
    return tuple(max(1, int(math.ceil(2 * math.pi * x / kspacing - 1e-9))) for x in b)


def parse_kpoints(text: str, lattice_ang: np.ndarray):
    """Return ('mesh', mesh, rmg_shift_flags, description) or ('list', kfrac, weights, desc)."""
    lines = [l.strip() for l in text.splitlines()]
    if len(lines) < 3:
        raise ConversionError("KPOINTS is too short")
    try:
        n = int(lines[1].split()[0])
    except (ValueError, IndexError):
        raise ConversionError("bad KPOINTS second line %r" % lines[1])
    mode = lines[2][:1].lower()
    if n == 0:
        if mode == "a":
            rk = float(lines[3].split()[0])
            mesh = mesh_from_length(lattice_ang, rk)
            return ("mesh", mesh, (0, 0, 0), "fully automatic, R_k = %g" % rk)
        if mode in ("g", "m"):
            mesh = tuple(int(x) for x in lines[3].split()[:3])
            if len(mesh) != 3 or min(mesh) < 1:
                raise ConversionError("bad KPOINTS mesh line %r" % lines[3])
            shift = (0.0, 0.0, 0.0)
            if len(lines) > 4 and lines[4]:
                shift = tuple(float(x) for x in lines[4].split()[:3])
            flags = _shift_flags(mesh, shift, mode == "m")
            return ("mesh", mesh, flags, ("Gamma-centred" if mode == "g" else "Monkhorst-Pack"))
        raise ConversionError("KPOINTS automatic mode with generating vectors is not supported")
    if mode == "l":
        raise ConversionError("line-mode KPOINTS (band paths) is not converted; use an RMG "
                              "kpoints_bandstructure block")
    if mode in ("c", "k"):
        raise ConversionError("Cartesian explicit KPOINTS lists are not converted; give them in "
                              "reciprocal coordinates")
    rows = [[float(x) for x in l.split()[:4]] for l in lines[3:3 + n]]
    if len(rows) != n or any(len(r) < 4 for r in rows):
        raise ConversionError("explicit KPOINTS list is incomplete")
    return ("list", np.array([r[:3] for r in rows]), [r[3] for r in rows], "explicit list")


@dataclass
class PotcarInfo:
    titles: List[str]
    elements: List[str]
    zval: List[float]
    enmax: List[float]


def parse_potcar(text: str) -> PotcarInfo:
    """Species order, ZVAL and ENMAX of a (concatenated) POTCAR. Nothing else is used."""
    titles = re.findall(r"TITEL\s*=\s*(.+)", text)
    zval = [float(x) for x in re.findall(r"ZVAL\s*=\s*([-+0-9.Ee]+)", text)]
    enmax = [float(x) for x in re.findall(r"ENMAX\s*=\s*([-+0-9.Ee]+)", text)]
    if not titles or len(titles) != len(zval):
        raise ConversionError("POTCAR: could not find one TITEL and ZVAL per species")
    elements = []
    for t in titles:
        f = t.split()
        el = element_from_label(f[1] if len(f) > 1 else f[0])
        if el is None:
            raise ConversionError("POTCAR: cannot tell the element of %r" % t)
        elements.append(el)
    return PotcarInfo([t.strip() for t in titles], elements, zval, enmax)


# ---------------------------------------------------------------------------------------------
# Tag tables
# ---------------------------------------------------------------------------------------------

GGA_MAP = {"PE": "pbe", "PS": "pbesol", "91": "pw91", "RE": "revpbe", "CA": "pz", "PZ": "pz"}
METAGGA_MAP = {"SCAN": "scan", "R2SCAN": "r2scan", "TPSS": "tpss", "M06L": "m06l", "MBJ": "tb09"}
LIBXC_ONLY = {"scan", "r2scan"}

HANDLED = {"ENCUT", "ISPIN", "MAGMOM", "NUPDOWN", "LHFCALC", "HFSCREEN", "AEXX", "GGA",
           "METAGGA", "ISMEAR", "SIGMA", "NBANDS", "NELECT", "EDIFF", "EDIFFG", "NELM", "IBRION",
           "NSW", "ISIF", "SMASS", "TEBEG", "TEEND", "POTIM", "MDALGO", "ISYM", "AMIX", "ICHARG",
           "KSPACING", "KGAMMA", "NKRED", "NKREDX", "NKREDY", "NKREDZ", "IVDW", "LNONCOLLINEAR",
           "LSORBIT", "SAXIS", "SYSTEM", "LDAU", "ALDAC", "AGGAX", "AGGAC", "AMGGAX", "AMGGAC",
           "LUSE_VDW", "PRECFOCK", "GGA_COMPAT"}
RUN_CONTROL = {"PREC", "ALGO", "LREAL", "LWAVE", "LCHARG", "LVTOT", "LVHAR", "LORBIT", "NPAR",
               "NCORE", "KPAR", "NSIM", "LPLANE", "LSCALAPACK", "LSCALU", "ISTART", "NWRITE",
               "ADDGRID", "ENAUG", "NELMIN", "NELMDL", "BMIX", "AMIX_MAG", "BMIX_MAG", "IMIX",
               "MAXMIX", "LMAXMIX", "NBLOCK", "KBLOCK", "ROPT", "EMIN", "EMAX", "NEDOS", "LELF",
               "LAECHG", "LMIXTAU", "LSUBROT", "TIME", "WEIMIN", "EBREAK", "DEPER", "LMAXFOCK",
               "NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF", "LASPH", "IWAVPR", "ISTART",
               "LREAL_COMPAT", "NELMALL", "LPARD", "NBMOD", "EINT", "LSEPB", "LSEPK", "IOPT",
               "LLCALL", "NFREE", "ENCUTFOCK", "LLRAUG", "NOMEGA", "KPOINT_BSE"}
PHYSICS_UNSUPPORTED = {"LEPSILON", "LOPTICS", "LCALCEPS", "EFIELD", "EFIELD_PEAD", "LDIPOL",
                       "IDIPOL", "LCALCPOL", "LBERRY", "LSOL", "EB_K", "LCHIMAG", "ICORELEVEL",
                       "IMAGES", "LAMBDA", "I_CONSTRAINED_M", "M_CONSTR", "LRPA", "LHYPERFINE",
                       "LNMR_SYM_RED", "IALGO", "HFRCUT", "LMODELHF", "LTHOMAS", "LSCAAWARE",
                       "LVDW", "VDW_S6", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ", "LORBMOM",
                       "QSPIRAL", "LSPIRAL", "LZEROZ", "LAMBDA_D_K", "LUSE_VDW"}


# ---------------------------------------------------------------------------------------------
# Converter
# ---------------------------------------------------------------------------------------------


def _review_tags(incar: Dict[str, str], notes: Notes) -> None:
    ignored = []
    for k in incar:
        if k in PHYSICS_UNSUPPORTED:
            notes.unsupported("INCAR %s is not converted" % k)
        elif k in HANDLED:
            continue
        elif k in RUN_CONTROL:
            ignored.append(k)
        else:
            notes.warn("INCAR %s is not recognised by the converter and was ignored" % k)
    if ignored:
        notes.info("run-control INCAR tags with no RMG counterpart ignored: " + " ".join(ignored))


def convert(poscar: str, incar: Optional[str] = None, kpoints: Optional[str] = None,
            species: Optional[Sequence[str]] = None, potcar: Optional[str] = None,
            pseudos: Optional[Dict[str, str]] = None, system_charge: Optional[float] = None,
            grid: Optional[Sequence[int]] = None, grid_spacing: Optional[float] = None,
            grid_multiple: int = 4, lenient: bool = False) -> Tuple[RmgDeck, Notes]:
    """Convert VASP input texts. ``pseudos`` maps element -> UPF path."""
    notes = Notes(lenient)
    pinfo = parse_potcar(potcar) if potcar else None
    pos = parse_poscar(poscar, species, pinfo.elements if pinfo is not None else None)
    tags = parse_incar(incar) if incar else {}
    _review_tags(tags, notes)
    elems_order = []
    for s in pos.symbols:
        if s not in elems_order:
            elems_order.append(s)
    if pinfo is not None and pinfo.elements != elems_order:
        raise ConversionError("POTCAR species %s do not match the POSCAR species %s"
                              % (pinfo.elements, elems_order))

    deck = RmgDeck()
    deck.set("description", sanitize_description("converted from VASP input " + pos.comment))

    def has(k):
        return k in tags

    def fval(k, d=None):
        return incar_float(tags[k]) if k in tags else d

    def ival(k, d=None):
        return incar_int(tags[k]) if k in tags else d

    def bval(k, d=False):
        return incar_bool(tags[k]) if k in tags else d

    # -- calculation mode ------------------------------------------------------------------
    nsw = ival("NSW", 0)
    ibrion = ival("IBRION", -1 if nsw in (0, 1) else 0)
    isif = ival("ISIF", 0 if ibrion == 0 else 2)
    if ival("ICHARG", 2) == 11:
        mode = "NSCF"
    elif nsw <= 0 or ibrion == -1:
        mode = "Quench Electrons"
        if nsw > 0:
            notes.info("IBRION = -1 with NSW > 0 is a static calculation in RMG")
    elif ibrion in (1, 2, 3):
        mode = "Relax Structure"
        deck.set("max_md_steps", nsw)
        if isif == 3:
            deck.set("cell_relax", True)
        elif isif > 3:
            notes.unsupported("ISIF = %d (constrained cell relaxation) is not converted" % isif)
        ediffg = fval("EDIFFG", 10 * fval("EDIFF", 1e-4))
        if ediffg < 0:
            deck.set("relax_max_force", abs(ediffg) / HA_EV / BOHR_PER_ANGSTROM)
        else:
            notes.warn("EDIFFG > 0 (energy criterion) is not converted; RMG stops relaxations on "
                       "relax_max_force (default kept)")
        notes.info("IBRION = %d: RMG's default relax_method (LBFGS) is used" % ibrion)
    elif ibrion == 0:
        mdalgo = ival("MDALGO", 0)
        smass = fval("SMASS", -3.0)
        if mdalgo not in (0, 2):
            notes.unsupported("MDALGO = %d is not converted" % mdalgo)
        if smass == -3 and mdalgo == 0:
            mode = "Constant Volume And Energy"
        elif smass >= 0 or mdalgo == 2:
            mode = "Constant Temperature And Energy"
            deck.set("md_temperature_control", "Nose Hoover Chains")
            deck.set("md_temperature", fval("TEBEG", 0.0))
            if has("TEEND") and fval("TEEND") != fval("TEBEG", 0.0):
                notes.unsupported("a temperature ramp (TEEND != TEBEG) is not converted")
        else:
            notes.unsupported("SMASS = %g is not converted" % smass)
            mode = "Constant Volume And Energy"
        deck.set("max_md_steps", nsw)
        notes.warn("POTIM is not converted; set ionic_time_step for RMG by hand")
    else:
        notes.unsupported("IBRION = %d is not converted" % ibrion)
        mode = "Quench Electrons"
    deck.set("calculation_mode", mode)

    # -- pseudopotentials ------------------------------------------------------------------
    zval_upf: Dict[str, float] = {}
    pp_dir = None
    if pseudos:
        missing = [e for e in elems_order if e not in pseudos]
        if missing:
            raise ConversionError("--pseudo given for some elements but not for %s" % missing)
        dirs = {os.path.dirname(pseudos[e]) for e in elems_order}
        if len(dirs) != 1:
            raise ConversionError("all --pseudo files must be in one directory (RMG has a single "
                                  "pseudo_dir)")
        for e in elems_order:
            p = pseudos[e]
            if os.path.isfile(p):
                h = read_upf_header(p)
                check_upf_for_rmg(h, notes)
                if h.element and element_from_label(h.element) != e:
                    raise ConversionError("%s is a pseudopotential for %s, not %s"
                                          % (p, h.element, e))
                if h.z_valence is not None:
                    zval_upf[e] = h.z_valence
            else:
                notes.warn("pseudopotential %s not found; its valence charge and format could "
                           "not be checked" % p)
        pp_dir = dirs.pop() or "."
    else:
        notes.warn("POTCAR (PAW) cannot be converted: no pseudopotential block is written, so "
                   "RMG uses its internal norm-conserving pseudopotentials "
                   "(internal_pseudo_type default). Pass --pseudo El=file.upf to choose UPF files")
    nel_known = all(e in zval_upf for e in elems_order)
    zsum = sum(zval_upf[s] for s in pos.symbols) if nel_known else None

    # -- charge ----------------------------------------------------------------------------
    charge = 0.0
    if system_charge is not None:
        charge = float(system_charge)
        if has("NELECT"):
            notes.info("--system-charge overrides NELECT")
    elif has("NELECT"):
        if pinfo is None:
            notes.unsupported("NELECT is relative to the POTCAR valences, which are unknown "
                              "without a POTCAR; give --potcar or --system-charge")
        else:
            zpot = dict(zip(pinfo.elements, pinfo.zval))
            charge = sum(zpot[s] for s in pos.symbols) - fval("NELECT")
            notes.info("system_charge = sum(POTCAR ZVAL) - NELECT = %g" % charge)
    if abs(charge) > 1e-12:
        deck.set("system_charge", charge)
    nel = None if zsum is None else zsum - charge

    # -- functional ------------------------------------------------------------------------
    lhf = bval("LHFCALC")
    gga = tags.get("GGA", "").strip().upper()
    mgga = tags.get("METAGGA", "").strip().upper()
    xc = None
    hybrid = False
    if mgga and mgga not in ("NONE", "--"):
        if lhf:
            notes.unsupported("meta-GGA hybrids are not converted")
        elif mgga not in METAGGA_MAP:
            notes.unsupported("METAGGA = %s has no RMG counterpart" % mgga)
        else:
            xc = METAGGA_MAP[mgga]
    elif lhf:
        hybrid = True
        custom = [k for k in ("ALDAC", "AGGAX", "AGGAC", "AMGGAX", "AMGGAC") if has(k)]
        if gga == "B3":
            xc = "b3lyp"
            if has("AEXX") and abs(fval("AEXX") - 0.2) > 1e-9:
                deck.set("exx_fraction", fval("AEXX"))
                notes.warn("B3LYP with AEXX = %g: only the exact-exchange fraction is changed"
                           % fval("AEXX"))
        elif gga in ("", "PE"):
            if custom:
                notes.unsupported("custom hybrid mixing (%s) is not converted" % " ".join(custom))
            hfs = fval("HFSCREEN", 0.0)
            if hfs == 0.0:
                xc = "pbe0"
            elif 0.19 <= hfs <= 0.21:
                xc = "hse"
                notes.info("HFSCREEN = %g A^-1 (%.4f bohr^-1) mapped to RMG HSE06, which uses "
                           "0.106 bohr^-1" % (hfs, hfs / BOHR_PER_ANGSTROM))
            else:
                notes.unsupported("HFSCREEN = %g A^-1 is not HSE06; RMG has no input for the "
                                  "screening length" % hfs)
            if has("AEXX") and abs(fval("AEXX") - 0.25) > 1e-9:
                deck.set("exx_fraction", fval("AEXX"))
        else:
            notes.unsupported("LHFCALC with GGA = %s is not converted" % gga)
    elif gga:
        if gga not in GGA_MAP:
            notes.unsupported("GGA = %s has no RMG counterpart in this converter" % gga)
        else:
            xc = GGA_MAP[gga]
    else:
        notes.info("GGA not set: VASP takes the functional from the POTCAR, RMG from the "
                   "pseudopotentials (AUTO_XC)")
    if xc:
        deck.set("exchange_correlation_type", xc)
        if xc in LIBXC_ONLY:
            notes.warn("exchange_correlation_type %s needs an RMG build with libxc" % xc)
    if has("GGA_COMPAT"):
        notes.info("GGA_COMPAT not converted")
    if has("LDAU") and bval("LDAU"):
        notes.unsupported("DFT+U (LDAU) is not converted")

    # -- spin ------------------------------------------------------------------------------
    ispin = ival("ISPIN", 1)
    noncoll = bval("LNONCOLLINEAR") or bval("LSORBIT")
    nat = len(pos.symbols)
    spins = [()] * nat
    spin_cols = 0

    def ratio_of(moment: float, el: str) -> float:
        if el in zval_upf:
            return moment / (2.0 * zval_upf[el])
        return moment / 10.0

    if noncoll:
        if has("SAXIS") and [incar_float(x) for x in tags["SAXIS"].split()[:3]] != [0.0, 0.0, 1.0]:
            notes.unsupported("SAXIS other than 0 0 1 is not converted")
        mag = expand_magmom(tags["MAGMOM"]) if has("MAGMOM") else [1.0] * (3 * nat)
        if len(mag) != 3 * nat:
            raise ConversionError("noncollinear MAGMOM needs 3 x %d values, got %d"
                                  % (nat, len(mag)))
        spins = []
        for i, el in enumerate(pos.symbols):
            mx, my, mz = mag[3 * i:3 * i + 3]
            m = math.sqrt(mx * mx + my * my + mz * mz)
            th = math.degrees(math.acos(mz / m)) if m > 0 else 0.0
            ph = math.degrees(math.atan2(my, mx)) if m > 0 else 0.0
            spins.append((clamp_spin_ratio(ratio_of(m, el), notes, "atom %d" % (i + 1)), th, ph))
        spin_cols = 3
        deck.set("noncollinear", True)
        if bval("LSORBIT"):
            deck.set("spinorbit", True)
            notes.warn("spin-orbit needs fully relativistic pseudopotentials in RMG")
    elif ispin == 2:
        mag = expand_magmom(tags["MAGMOM"]) if has("MAGMOM") else [1.0] * nat
        if len(mag) != nat:
            raise ConversionError("MAGMOM has %d values for %d atoms" % (len(mag), nat))
        spins = [(clamp_spin_ratio(ratio_of(m, el), notes, "atom %d" % (i + 1)),)
                 for i, (m, el) in enumerate(zip(mag, pos.symbols))]
        spin_cols = 1
        deck.set("spin_polarization", True)
    elif ispin != 1:
        raise ConversionError("ISPIN must be 1 or 2")
    if spin_cols:
        if nel_known:
            notes.info("RMG initial spin column = MAGMOM / (2 * UPF z_valence)")
        else:
            notes.warn("RMG initial spin column = MAGMOM / 10 (clamped to +-0.5): the valence "
                       "charges are unknown without --pseudo; it is only the starting guess")
    if has("NUPDOWN") and fval("NUPDOWN") >= 0:
        if ispin != 2:
            raise ConversionError("NUPDOWN requires ISPIN = 2")
        deck.set("tot_magnetization", fval("NUPDOWN"))
        notes.info("tot_magnetization is a recent RMG key; older RMG builds ignore it silently")

    # -- occupations -----------------------------------------------------------------------
    ismear = ival("ISMEAR", 1)
    sigma = fval("SIGMA", 0.2)
    if ismear == -5:
        deck.set("occupations_type", "Tetrahedron")
        deck.set("tetra_method", "Bloechl")
    elif ismear == -4:
        deck.set("occupations_type", "Tetrahedron")
        deck.set("tetra_method", "Linear")
    elif ismear in (-1, 0) or 1 <= ismear <= 5:
        if sigma > 2.0:
            notes.unsupported("SIGMA = %g eV exceeds RMG's 2 eV limit" % sigma)
        if ismear == -1:
            deck.set("occupations_type", "Fermi Dirac")
        else:
            deck.set("occupations_type", "MethfesselPaxton")
            deck.set("MP_order", max(ismear, 0))
            if ismear == 0:
                notes.info("VASP Gaussian smearing = RMG MethfesselPaxton with MP_order 0 "
                           "(RMG's 'Gaussian' occupation is a different function)")
        deck.set("occupation_electron_temperature_eV", sigma)
    else:
        notes.unsupported("ISMEAR = %d is not converted" % ismear)

    if has("NBANDS"):
        if nel is None:
            notes.warn("NBANDS not converted: it needs the valence charge of each RMG "
                       "pseudopotential (give --pseudo); RMG's default is used")
        else:
            nocc = nocc_rmg(nel, noncoll)
            unocc = ival("NBANDS") - nocc
            if unocc < 0:
                notes.warn("NBANDS = %d is below RMG's %d occupied states; set to 0 unoccupied"
                           % (ival("NBANDS"), nocc))
                unocc = 0
            deck.set("unoccupied_states_per_kpoint", unocc)

    # -- k-points --------------------------------------------------------------------------
    if kpoints:
        kp = parse_kpoints(kpoints, pos.lattice)
    else:
        ksp = fval("KSPACING", 0.5)
        mesh = mesh_from_kspacing(pos.lattice, ksp)
        flags = (0, 0, 0) if bval("KGAMMA", True) else _shift_flags(mesh, (0, 0, 0), True)
        kp = ("mesh", mesh, flags, "from KSPACING = %g" % ksp)
        if not has("KSPACING"):
            notes.warn("no KPOINTS file and no KSPACING: VASP's default KSPACING = 0.5 used")
    notes.info("k-points: %s" % kp[3])

    lattice_bohr = pos.lattice * BOHR_PER_ANGSTROM
    st = Structure(lattice_bohr, pos.symbols, pos.frac, pos.selective, spins)
    st, orient = apply_orientation(st, notes)
    if kp[0] == "mesh":
        kmesh = kp[1]
        deck.set("kpoint_mesh", list(kmesh))
        deck.set("kpoint_is_shift", list(kp[2]))
    else:
        kmesh = None
        write_kpoint_list(deck, np.asarray(kp[1]) @ orient.basis_change.T, kp[2])

    # -- exact exchange --------------------------------------------------------------------
    if hybrid:
        if kmesh is None:
            notes.unsupported("hybrid functionals with an explicit k-point list are not converted")
        else:
            nk = ival("NKRED", 1)
            red = [ival("NKREDX", nk), ival("NKREDY", nk), ival("NKREDZ", nk)]
            for k, r in zip(kmesh, red):
                if k % r:
                    raise ConversionError("NKRED %s does not divide the k mesh %s"
                                          % (red, list(kmesh)))
            q = [k // r for k, r in zip(kmesh, red)]
            check_exx_qmesh(kmesh, q)
            deck.set("exx_qpoint_mesh", q)
            notes.info("exx_qpoint_mesh = k mesh / NKRED; it needs an RMG build with "
                       "exx_qpoint_mesh support")
        notes.warn("VASP's treatment of the Coulomb singularity is not RMG's; RMG defaults "
                   "(exxdiv_treatment gygi-baldereschi, x_gamma_extrapolation true) are used")
        if has("PRECFOCK"):
            notes.info("PRECFOCK not converted")

    # -- vdW -------------------------------------------------------------------------------
    ivdw = ival("IVDW", 0)
    if ivdw in (1, 10):
        deck.set("vdw_corr", "DFT-D2")
    elif ivdw == 11:
        deck.set("vdw_corr", "DFT-D3")
        deck.set("dftd3_version", 3)
    elif ivdw == 12:
        deck.set("vdw_corr", "DFT-D3")
        deck.set("dftd3_version", 4)
    elif ivdw != 0:
        notes.unsupported("IVDW = %d is not converted" % ivdw)

    # -- symmetry --------------------------------------------------------------------------
    isym = ival("ISYM", None)
    if isym is not None and isym <= 0:
        deck.set("use_symmetry", 0)
        if isym == -1:
            deck.set("time_reversal", False)

    # -- grid ------------------------------------------------------------------------------
    ecut_ev = fval("ENCUT")
    if ecut_ev is None and pinfo is not None and pinfo.enmax:
        ecut_ev = max(pinfo.enmax)
        notes.info("ENCUT not set: VASP's default max(ENMAX) = %g eV from the POTCAR used" % ecut_ev)
    if ecut_ev is None and grid is None and grid_spacing is None:
        raise ConversionError("ENCUT is not set and there is no POTCAR: give --grid or "
                              "--grid-spacing")
    g = wavefunction_grid(st.lattice, ecut_ev / RY_EV if ecut_ev else None, grid_spacing, grid,
                          grid_multiple, notes)
    deck.set("wavefunction_grid", list(g.grid))
    notes.info("wavefunction_grid %s %s; spacings %.4f %.4f %.4f bohr; equivalent cutoff "
               "%.2f Ry" % (" ".join(map(str, g.grid)), g.how.replace("ecutwfc", "ENCUT"),
                            *g.spacing_bohr, g.ecut_equiv_ry))
    if grid is None and grid_spacing is None:
        notes.warn("ENCUT is a PAW cutoff; RMG's norm-conserving pseudopotentials usually need "
                   "a finer grid. Check convergence or use --grid-spacing")

    # -- SCF control -----------------------------------------------------------------------
    set_scf_thresholds(deck, fval("EDIFF", 1e-4) / HA_EV, notes, hybrid=hybrid, source="EDIFF",
                       spin=(ispin == 2 or noncoll))
    deck.set("max_scf_steps", ival("NELM", 60))
    if has("AMIX"):
        deck.set("charge_density_mixing", fval("AMIX"))

    if pp_dir is not None:
        deck.set("pseudo_dir", pp_dir)
        deck.set_block("pseudopotential",
                       ["  %s  %s" % (e, os.path.basename(pseudos[e])) for e in elems_order])
    write_structure(deck, st, spin_cols)
    return deck, notes


def _read(path: Optional[str]) -> Optional[str]:
    if path is None:
        return None
    with open(path) as f:
        return f.read()


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog=PROG, description=__doc__.split("\n\n")[0])
    ap.add_argument("-d", "--dir", default=".", help="directory with the VASP files")
    ap.add_argument("--poscar", help="POSCAR/CONTCAR (default DIR/POSCAR, else DIR/CONTCAR)")
    ap.add_argument("--incar", help="INCAR (default DIR/INCAR if present)")
    ap.add_argument("--kpoints", help="KPOINTS (default DIR/KPOINTS if present)")
    ap.add_argument("--potcar", help="POTCAR, read for TITEL/ZVAL/ENMAX only (default "
                                     "DIR/POTCAR if present)")
    ap.add_argument("--species", help="species names for a VASP 4 POSCAR, e.g. \"Fe Si\"")
    ap.add_argument("--pseudo", action="append", default=[], metavar="El=FILE",
                    help="UPF file for an element (repeat for every element)")
    ap.add_argument("--system-charge", type=float, help="net charge (overrides NELECT)")
    ap.add_argument("-o", "--output", default="-", help="RMG input file (default stdout)")
    ap.add_argument("--notes", default=None,
                    help="notes file (default <output>.notes; none when writing to stdout)")
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument("--grid", nargs=3, type=int, metavar=("NX", "NY", "NZ"))
    grp.add_argument("--grid-spacing", type=float, metavar="H", help="target spacing in bohr")
    ap.add_argument("--grid-multiple", type=int, default=4)
    ap.add_argument("--lenient", action="store_true",
                    help="drop unsupported features with a warning instead of stopping")
    a = ap.parse_args(argv)

    def pick(explicit, *names):
        if explicit:
            return explicit
        for n in names:
            p = os.path.join(a.dir, n)
            if os.path.isfile(p):
                return p
        return None

    try:
        poscar = pick(a.poscar, "POSCAR", "CONTCAR")
        if poscar is None:
            raise ConversionError("no POSCAR or CONTCAR found")
        pseudos = {}
        for item in a.pseudo:
            if "=" not in item:
                raise ConversionError("--pseudo needs El=FILE, got %r" % item)
            el, path = item.split("=", 1)
            pseudos[el.strip()] = path.strip()
        deck, notes = convert(_read(poscar), _read(pick(a.incar, "INCAR")),
                              _read(pick(a.kpoints, "KPOINTS")),
                              a.species.split() if a.species else None,
                              _read(pick(a.potcar, "POTCAR")), pseudos or None,
                              a.system_charge, a.grid, a.grid_spacing, a.grid_multiple,
                              a.lenient)
        write_outputs(PROG, deck, notes, a.output, a.notes,
                      "%s notes for %s (RMG decks cannot hold comments)" % (PROG, poscar))
    except ConversionError as exc:
        print("%s: error: %s" % (PROG, exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
