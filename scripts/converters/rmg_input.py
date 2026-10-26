"""Shared building blocks for the QE -> RMG and VASP -> RMG input converters.

This module holds everything that does not depend on the source code's input format:

* :class:`RmgDeck`, an ordered model of an RMG input file and its writer, with the same
  restrictions RMG's own parser imposes (``Input/RmgInputFile.cpp``): every entry is
  ``key = "value"``, a value may span several lines inside the quotes, a value may not
  contain ``=`` or ``"``, and ``#`` starts a comment, so it may not appear in a value either.
  The converters never write comment lines; their notes go to stderr and a ``.notes`` file.
* :func:`parse_rmg_deck` and :func:`deck_structure`, a small reader used by the tests to read
  a written deck back and compare cell and positions numerically.
* Unit constants taken from RMG's ``Headers/const.h`` so that a deck written in bohr means
  exactly what RMG would compute from the same input in Angstrom.
* The real-space grid mapping from a plane-wave cutoff (:func:`wavefunction_grid`).
* A Python copy of RMG's lattice re-orientation (``Lattice::lat2ibrav`` and
  ``Lattice::rotate_vectors`` in ``RmgLib/src/Lattice.cpp``) and :func:`orient_for_rmg`,
  which rewrites a cell and its fractional coordinates so that the structure RMG builds is
  the structure that was asked for.
* A minimal UPF header reader (valence charge, element, pseudopotential type).

Only the Python standard library and numpy are used.
"""

from __future__ import annotations

import itertools
import math
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------------------------
# Units (RMG Headers/const.h)
# ---------------------------------------------------------------------------------------------

#: bohr per Angstrom, RMG's ``A_a0``.
BOHR_PER_ANGSTROM = 1.88972612499359
#: Hartree in eV, RMG's ``Ha_eV``.
HA_EV = 27.2113845
#: Rydberg in eV.
RY_EV = HA_EV / 2.0


def angstrom_to_bohr(x):
    """Convert a length (or array of lengths) from Angstrom to bohr with RMG's constant."""
    return np.asarray(x, dtype=float) * BOHR_PER_ANGSTROM


def ev_to_ry(x: float) -> float:
    return float(x) / RY_EV


def ry_to_ev(x: float) -> float:
    return float(x) * RY_EV


# ---------------------------------------------------------------------------------------------
# Errors and notes
# ---------------------------------------------------------------------------------------------


class ConversionError(Exception):
    """The input asks for something the converter cannot represent faithfully in RMG."""


class Notes:
    """Collects what a conversion changed, dropped or could not check.

    ``warn`` is for anything the user must look at (a default that differs, a value that was
    clamped, a quantity that could not be derived). ``info`` is for plain bookkeeping.
    ``unsupported`` refuses the conversion with :class:`ConversionError` unless the converter
    runs in lenient mode, in which case the dropped feature is recorded as a warning.
    """

    def __init__(self, lenient: bool = False):
        self.lenient = lenient
        self.warnings: List[str] = []
        self.infos: List[str] = []

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)

    def info(self, msg: str) -> None:
        self.infos.append(msg)

    def unsupported(self, msg: str) -> None:
        if not self.lenient:
            raise ConversionError(msg + " (use --lenient to drop it and convert anyway)")
        self.warnings.append("DROPPED: " + msg)

    def text(self, header: str = "") -> str:
        lines = []
        if header:
            lines.append(header)
        for w in self.warnings:
            lines.append("WARNING: " + w)
        for i in self.infos:
            lines.append("NOTE: " + i)
        return "\n".join(lines) + ("\n" if lines else "")

    def emit(self, prog: str, stream=None) -> None:
        stream = sys.stderr if stream is None else stream
        for w in self.warnings:
            print("%s: warning: %s" % (prog, w), file=stream)
        for i in self.infos:
            print("%s: note: %s" % (prog, i), file=stream)


# ---------------------------------------------------------------------------------------------
# The RMG deck
# ---------------------------------------------------------------------------------------------

_KEY_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_FORBIDDEN_IN_VALUE = ("=", '"', "#")


def fmt_float(x: float, digits: int = 12) -> str:
    """Shortest exact-enough text for a float (``%.<digits>g``, no trailing noise)."""
    x = float(x)
    if x == 0.0:
        return "0.0"
    s = "%.*g" % (digits, x)
    if "e" not in s and "." not in s and "inf" not in s and "nan" not in s:
        s += ".0"
    return s


def fmt_bool(b: bool) -> str:
    return "true" if b else "false"


def fmt_ints(v: Iterable[int]) -> str:
    return " ".join(str(int(i)) for i in v)


class RmgDeck:
    """Ordered ``key -> value`` model of an RMG input file.

    Values are stored as text. Keys and values are validated when they are set, so an
    invalid deck can never be rendered.
    """

    def __init__(self):
        self._entries: Dict[str, str] = {}

    # -- validation ----------------------------------------------------------------------
    @staticmethod
    def validate_key(key: str) -> None:
        if not _KEY_RE.match(key):
            raise ValueError("invalid RMG key %r" % key)

    @staticmethod
    def validate_value(key: str, value: str) -> None:
        for ch in _FORBIDDEN_IN_VALUE:
            if ch in value:
                raise ValueError(
                    "value of RMG key %r contains %r, which RMG's input parser rejects or "
                    "treats as a comment: %r" % (key, ch, value))
        for line in value.splitlines():
            if line.strip().startswith("#"):
                raise ValueError("value of RMG key %r has a line starting with '#'" % key)

    # -- setters -------------------------------------------------------------------------
    def set(self, key: str, value) -> None:
        """Set a single-line value. Booleans, ints, floats and int triplets are formatted."""
        if isinstance(value, bool):
            text = fmt_bool(value)
        elif isinstance(value, (int, np.integer)):
            text = str(int(value))
        elif isinstance(value, (float, np.floating)):
            text = fmt_float(value)
        elif isinstance(value, (list, tuple, np.ndarray)):
            text = " ".join(fmt_float(v) if isinstance(v, (float, np.floating)) else str(v)
                            for v in value)
        else:
            text = str(value)
        if "\n" in text:
            raise ValueError("use set_block for multi-line value of %r" % key)
        self.validate_key(key)
        self.validate_value(key, text)
        self._entries[key] = text

    def set_block(self, key: str, lines: Sequence[str]) -> None:
        """Set a multi-line value (lattice vectors, atoms, pseudopotentials, k-points)."""
        self.validate_key(key)
        text = "\n".join(str(l).rstrip() for l in lines)
        self.validate_value(key, text)
        self._entries[key] = "\n" + text + "\n"

    def get(self, key: str, default=None):
        v = self._entries.get(key, default)
        return v.strip() if isinstance(v, str) else v

    def __contains__(self, key: str) -> bool:
        return key in self._entries

    def remove(self, key: str) -> None:
        self._entries.pop(key, None)

    def keys(self) -> List[str]:
        return list(self._entries)

    # -- writer --------------------------------------------------------------------------
    def render(self) -> str:
        out = []
        for key, value in self._entries.items():
            out.append('%s = "%s"' % (key, value))
        text = "\n".join(out) + "\n"
        check_rmg_text(text)
        return text

    def write(self, path: str) -> None:
        with open(path, "w") as f:
            f.write(self.render())


def check_rmg_text(text: str) -> None:
    """Final guard on the rendered text: no comment lines and every entry parses back."""
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            raise ValueError("rendered deck contains a comment line: %r" % line)
    stripped = re.sub(r'"[^"]*"', '""', text, flags=re.S)
    for line in stripped.splitlines():
        if line.strip() and not re.match(r'^\s*[A-Za-z_]\w*\s*=\s*""\s*$', line):
            raise ValueError("rendered deck has a malformed line: %r" % line)


def sanitize_description(text: str, limit: int = 200) -> str:
    """Make free text safe for the ``description`` key."""
    text = " ".join(str(text).split())
    for ch in _FORBIDDEN_IN_VALUE:
        text = text.replace(ch, " ")
    text = " ".join(text.split())
    return text[:limit]


# ---------------------------------------------------------------------------------------------
# Reading a deck back (tests and sanity checks)
# ---------------------------------------------------------------------------------------------


def parse_rmg_deck(text: str) -> Dict[str, str]:
    """Parse RMG input text into ``{key: value}`` the way RMG's preprocessor does.

    Lines starting with ``#`` are dropped, anything after a ``#`` is a trailing comment, and
    entries are ``key = "value"`` with values possibly spanning several lines.
    """
    kept = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "#" in line:
            line = line[: line.index("#")].strip()
        kept.append(line)
    body = "\n".join(kept)
    out: Dict[str, str] = {}
    pos = 0
    pat = re.compile(r'\s*([A-Za-z_]\w*)\s*=\s*"([^"]*)"', re.S)
    while pos < len(body):
        m = pat.match(body, pos)
        if not m:
            rest = body[pos:].strip()
            if not rest:
                break
            raise ValueError("malformed RMG input near %r" % rest[:60])
        out[m.group(1)] = m.group(2).strip()
        pos = m.end()
    return out


@dataclass
class DeckStructure:
    lattice_bohr: np.ndarray            # rows are lattice vectors
    symbols: List[str]
    cart_bohr: np.ndarray               # Cartesian positions
    frac: np.ndarray                    # fractional positions
    movable: np.ndarray                 # (N, 3) ints
    moments: List[Tuple[float, ...]]    # extra columns after the movable flags


def deck_structure(deck: Dict[str, str]) -> DeckStructure:
    """Rebuild the cell and the atoms from a parsed deck (lattice_vector form only)."""
    units = deck.get("lattice_units", "Bohr")
    vals = [float(x) for x in deck["lattice_vector"].split()]
    if len(vals) != 9:
        raise ValueError("lattice_vector must hold 9 numbers")
    lat = np.array(vals).reshape(3, 3)
    if units == "Angstrom":
        lat = lat * BOHR_PER_ANGSTROM
    elif units == "Alat":
        lat = lat * float(deck["a_length"])
    crd_type = deck.get("atomic_coordinate_type", "Absolute")
    crd_units = deck.get("crds_units", "Bohr")
    symbols, pos, mov, mom = [], [], [], []
    for line in deck["atoms"].splitlines():
        f = line.split()
        if not f:
            continue
        symbols.append(f[0])
        pos.append([float(x) for x in f[1:4]])
        rest = f[4:]
        if len(rest) in (1, 2):
            m = [int(rest[0] != "0")] * 3
            extra = tuple(float(x) for x in rest[1:])
        elif len(rest) in (3, 4, 6):
            m = [int(x != "0") for x in rest[:3]]
            extra = tuple(float(x) for x in rest[3:])
        elif len(rest) == 0:
            m, extra = [1, 1, 1], ()
        else:
            raise ValueError("bad atoms line %r" % line)
        mov.append(m)
        mom.append(extra)
    pos = np.array(pos, dtype=float)
    if crd_type == "Cell Relative":
        frac = pos
        cart = pos @ lat
    else:
        cart = pos * (BOHR_PER_ANGSTROM if crd_units == "Angstrom" else 1.0)
        frac = cart @ np.linalg.inv(lat)
    return DeckStructure(lat, symbols, cart, frac, np.array(mov, dtype=int), mom)


# ---------------------------------------------------------------------------------------------
# Structures
# ---------------------------------------------------------------------------------------------


@dataclass
class Structure:
    """A periodic structure in RMG's terms.

    ``lattice`` rows are the lattice vectors in bohr; ``frac`` are fractional coordinates;
    ``movable`` holds per-atom, per-axis 0/1 flags; ``spin`` holds per-atom RMG initial spin
    columns: ``()`` for none, ``(ratio,)`` for collinear, ``(ratio, theta, phi)`` for
    noncollinear. ``ratio`` is RMG's (up - down)/(up + down)/2 of the atomic density, in
    [-0.5, 0.5].
    """

    lattice: np.ndarray
    symbols: List[str]
    frac: np.ndarray
    movable: np.ndarray = None
    spin: List[Tuple[float, ...]] = None

    def __post_init__(self):
        self.lattice = np.asarray(self.lattice, dtype=float).reshape(3, 3)
        self.frac = np.asarray(self.frac, dtype=float).reshape(-1, 3)
        n = len(self.symbols)
        if self.frac.shape[0] != n:
            raise ValueError("number of symbols and positions differ")
        if self.movable is None:
            self.movable = np.ones((n, 3), dtype=int)
        self.movable = np.asarray(self.movable, dtype=int).reshape(n, 3)
        if self.spin is None:
            self.spin = [()] * n
        if abs(np.linalg.det(self.lattice)) < 1e-10:
            raise ConversionError("the lattice vectors are linearly dependent (zero volume)")

    @property
    def cart(self) -> np.ndarray:
        return self.frac @ self.lattice


def write_structure(deck: RmgDeck, st: Structure, spin_columns: int = 0) -> None:
    """Write lattice and atoms. ``spin_columns`` is 0, 1 (collinear) or 3 (noncollinear).

    The lattice is written in bohr and atoms in RMG "Cell Relative" (fractional) coordinates.
    Per-axis movable flags are written only when some atom is constrained on some axis but not
    all, otherwise the single movable column is used.
    """
    deck.set("lattice_units", "Bohr")
    deck.set_block("lattice_vector",
                   ["  %.12f  %.12f  %.12f" % tuple(v) for v in st.lattice])
    deck.set("atomic_coordinate_type", "Cell Relative")
    per_axis = spin_columns == 3 or any(len(set(m)) > 1 for m in st.movable.tolist())
    lines = []
    for i, sym in enumerate(st.symbols):
        x, y, z = st.frac[i]
        cols = ["  %-3s %16.12f %16.12f %16.12f" % (sym, x, y, z)]
        m = st.movable[i]
        if per_axis:
            cols.append("  %d %d %d" % tuple(m))
        else:
            cols.append("  %d" % int(m[0]))
        if spin_columns == 1:
            r = st.spin[i][0] if st.spin[i] else 0.0
            cols.append("  %+.6f" % r)
        elif spin_columns == 3:
            r, th, ph = st.spin[i] if st.spin[i] else (0.0, 0.0, 0.0)
            cols.append("  %+.6f %.6f %.6f" % (r, th, ph))
        lines.append("".join(cols))
    deck.set_block("atoms", lines)


def clamp_spin_ratio(r: float, notes: Optional[Notes] = None, what: str = "") -> float:
    """Clamp an initial spin ratio to RMG's allowed [-0.5, 0.5]."""
    if r > 0.5 or r < -0.5:
        c = max(-0.5, min(0.5, r))
        if notes is not None:
            notes.warn("initial spin ratio %.4f%s clamped to %.1f (RMG allows -0.5..0.5)"
                       % (r, (" for " + what) if what else "", c))
        return c
    return r


# ---------------------------------------------------------------------------------------------
# Grids
# ---------------------------------------------------------------------------------------------


def is_smooth(n: int, primes: Sequence[int] = (2, 3, 5, 7)) -> bool:
    for p in primes:
        while n % p == 0 and n > 1:
            n //= p
    return n == 1


def next_grid_size(n: int, multiple: int = 4, primes: Sequence[int] = (2, 3, 5, 7)) -> int:
    """Smallest m >= n that is a multiple of ``multiple`` and has only factors in ``primes``.

    RMG's own automatic grid (``Input/SetupGrids.cpp``) rounds every dimension to a multiple
    of 4, which keeps the multigrid levels usable; the 2,3,5,7-smooth condition keeps the FFTs
    (density interpolation, Poisson solver, exact exchange) on fast sizes.
    """
    m = max(int(n), 1)
    while True:
        if m % multiple == 0 and is_smooth(m, primes):
            return m
        m += 1


def min_points_for_cutoff(length_bohr: float, ecut_ry: float) -> int:
    """Minimum grid points along a lattice vector for a plane-wave cutoff (Ry).

    A grid of N points along a vector of length L has spacing h = L/N and resolves plane
    waves up to |G| = pi/h, i.e. a kinetic-energy cutoff of (pi/h)^2 Ry. This is exactly
    the "Equivalent energy cutoff (psi)" RMG prints (``Misc/WriteHeader.cpp``).
    """
    return int(math.ceil(length_bohr * math.sqrt(ecut_ry) / math.pi - 1e-9))


def equivalent_cutoff_ry(lattice_bohr: np.ndarray, grid: Sequence[int]) -> float:
    """Plane-wave cutoff (Ry) of the sphere that fits in the grid: (pi / max_i h_i)^2."""
    lengths = np.linalg.norm(np.asarray(lattice_bohr, dtype=float), axis=1)
    hmax = max(l / n for l, n in zip(lengths, grid))
    return (math.pi / hmax) ** 2


def planar_anisotropy(lattice_bohr: np.ndarray, grid: Sequence[int]) -> float:
    """RMG's planar anisotropy (``Misc/GetPlanarAnisotropy.cpp``); RMG stops above 1.5."""
    a = np.asarray(lattice_bohr, dtype=float)
    n = list(grid)
    dens = []
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        area = np.linalg.norm(np.cross(a[i], a[j]))
        dens.append(n[i] * n[j] / area)
    return max(dens) / min(dens)


@dataclass
class GridChoice:
    grid: Tuple[int, int, int]
    spacing_bohr: Tuple[float, float, float]
    ecut_equiv_ry: float
    how: str


def wavefunction_grid(lattice_bohr: np.ndarray, ecut_ry: Optional[float] = None,
                      spacing_bohr: Optional[float] = None,
                      grid: Optional[Sequence[int]] = None, multiple: int = 4,
                      notes: Optional[Notes] = None) -> GridChoice:
    """Choose RMG's ``wavefunction_grid``.

    Priority: an explicit ``grid``; else a target ``spacing_bohr``; else ``ecut_ry`` mapped
    with h = pi / sqrt(ecut_ry). Each dimension is rounded up with :func:`next_grid_size`.
    If the result exceeds RMG's planar-anisotropy limit (1.5) the coarsest directions are
    refined until it does not.
    """
    lat = np.asarray(lattice_bohr, dtype=float)
    lengths = np.linalg.norm(lat, axis=1)
    if grid is not None:
        g = tuple(int(x) for x in grid)
        if len(g) != 3 or min(g) < 1:
            raise ConversionError("--grid needs three positive integers")
        how = "set explicitly (--grid)"
        if notes is not None and planar_anisotropy(lat, g) > 1.5:
            notes.warn("the requested grid %s has planar anisotropy %.3f > 1.5; RMG will stop "
                       "with 'Planar Anisotropy too large'" % (g, planar_anisotropy(lat, g)))
    else:
        if spacing_bohr is not None:
            if spacing_bohr <= 0:
                raise ConversionError("--grid-spacing must be positive")
            nmin = [int(math.ceil(l / spacing_bohr - 1e-9)) for l in lengths]
            how = "from a target spacing of %.4f bohr" % spacing_bohr
        elif ecut_ry is not None:
            if ecut_ry <= 0:
                raise ConversionError("the wavefunction cutoff must be positive")
            nmin = [min_points_for_cutoff(l, ecut_ry) for l in lengths]
            how = ("from ecutwfc = %.4g Ry: h <= pi/sqrt(ecut) = %.5f bohr"
                   % (ecut_ry, math.pi / math.sqrt(ecut_ry)))
        else:
            raise ConversionError("no cutoff, grid spacing or grid given")
        g = [next_grid_size(n, multiple) for n in nmin]
        for _ in range(64):
            if planar_anisotropy(lat, g) <= 1.5:
                break
            h = [l / n for l, n in zip(lengths, g)]
            k = int(np.argmax(h))
            g[k] = next_grid_size(g[k] + 1, multiple)
        else:  # pragma: no cover - pathological cells only
            raise ConversionError("could not find a grid within RMG's planar anisotropy limit")
        if tuple(g) != tuple(next_grid_size(n, multiple) for n in nmin) and notes is not None:
            notes.warn("grid refined beyond the cutoff to stay under RMG's planar anisotropy "
                       "limit of 1.5")
        g = tuple(g)
    spacing = tuple(float(l / n) for l, n in zip(lengths, g))
    return GridChoice(g, spacing, equivalent_cutoff_ry(lat, g), how)


def potential_grid_refinement(ecutrho_ry: float, ecutwfc_ry: float) -> int:
    """RMG ``potential_grid_refinement`` from the density/wavefunction cutoff ratio.

    A fine grid r times denser than the wavefunction grid carries a density cutoff of
    r^2 times the wavefunction cutoff, so r = ceil(sqrt(ecutrho/ecutwfc)), at least 2
    (RMG's default, equivalent to QE's norm-conserving dual of 4) and at most 4 (RMG's limit).
    """
    dual = ecutrho_ry / ecutwfc_ry
    r = int(math.ceil(math.sqrt(dual) - 1e-9))
    return max(2, min(4, r))


# ---------------------------------------------------------------------------------------------
# Lattices: generic helpers and RMG's re-orientation
# ---------------------------------------------------------------------------------------------


def reciprocal(lattice: np.ndarray) -> np.ndarray:
    """Reciprocal vectors (rows) WITHOUT the 2*pi factor: a_i . b_j = delta_ij."""
    return np.linalg.inv(np.asarray(lattice, dtype=float)).T


def _eqq(x: float, y: float) -> bool:
    return abs(x - y) < 1.0e-5


def rmg_lat2ibrav(a: np.ndarray) -> int:
    """Python transcription of RMG's ``Lattice::lat2ibrav`` (RmgLib/src/Lattice.cpp)."""
    a0, a1, a2 = (np.asarray(v, dtype=float) for v in a)
    la, lb, lc = np.linalg.norm(a0), np.linalg.norm(a1), np.linalg.norm(a2)
    if abs(la - lb) < 1e-5:
        lb = la
    if abs(la - lc) < 1e-5:
        lc = la
    cab = float(a0 @ a1) / la / lb
    cac = float(a0 @ a2) / la / lc
    cbc = float(a2 @ a1) / lc / lb
    if abs(cab) < 1e-5:
        cab = 0.0
    if abs(cac) < 1e-5:
        cac = 0.0
    if abs(cbc) < 1e-5:
        cbc = 0.0
    if abs(cab - cac) < 1e-5:
        cac = cab
    if abs(cab - cbc) < 1e-5:
        cbc = cab
    e, ne = _eqq, (lambda x, y: not _eqq(x, y))
    ibrav = 14
    if e(cab, 0) and e(cac, 0) and e(cbc, 0):
        return 1 if (e(la, lb) and e(la, lc)) else 8
    if e(la, lb) and e(la, lc):
        if e(cab, cac) and e(cab, cbc):
            if e(cab, 0.0):
                ibrav = 1
            elif e(cab, 0.5):
                ibrav = 2
            elif e(cab, -1.0 / 3.0):
                ibrav = -3
            else:
                ibrav = 14
        elif e(cab, cac) and ne(cab, cbc):
            if e(abs(a0[0]), abs(a0[1])) and e(abs(a1[0]), abs(a1[1])):
                ibrav = 7
            else:
                ibrav = 3
        elif e(cab, -cac) and e(cab, cbc) and e(cab, 1.0 / 3.0):
            ibrav = 3
        elif e(abs(a0[0]), abs(a1[0])):
            ibrav = 14
    elif e(la, lb) and ne(la, lc):
        if e(cab, 0) and e(cac, 0) and e(cbc, 0):
            ibrav = 6
        elif e(cab, -0.5) and e(cac, 0) and e(cbc, 0):
            ibrav = 4
        elif e(cab, 0.5) and e(cac, 0) and e(cbc, 0):
            ibrav = 15
        elif e(cac, 0) and e(cbc, 0):
            ibrav = 12
        elif e(cac, -cbc):
            ibrav = 14
    elif e(la, lc) and ne(la, lb):
        if ne(cab, 0) and e(cac, 0) and e(cbc, 0):
            ibrav = 12
    elif e(lb, lc) and ne(la, lb):
        ibrav = 8 if (e(cab, 0) and e(cac, 0) and e(cbc, 0)) else 14
    else:
        if e(cab, 0) and e(cac, 0) and e(cbc, 0):
            ibrav = 8
        elif ne(cab, 0) and e(cac, 0) and e(cbc, 0):
            ibrav = 12
        else:
            ibrav = 14
    return ibrav


def rmg_rotate(a: np.ndarray) -> np.ndarray:
    """Python transcription of RMG's ``Lattice::rotate_vectors``.

    Returns the lattice vectors RMG will actually use when it is given ``a`` through
    ``lattice_vector`` (``Input/ReadCommon.cpp`` calls ``lat2ibrav`` then ``rotate_vectors``
    then keeps the result).
    """
    a = np.asarray(a, dtype=float)
    ib = rmg_lat2ibrav(a)
    if ib not in (2, 3, -3, 8, 1, 6, 7, 4, 15):
        return a.copy()
    m = np.linalg.norm(a, axis=1)
    r = np.zeros((3, 3))
    if ib in (1, 8, 6):
        r[0, 0], r[1, 1], r[2, 2] = m
    elif ib == 2:
        t = math.sqrt(2.0) * m[0] / 2.0
        r[0] = (t, t, 0)
        r[1] = (0, t, t)
        r[2] = (t, 0, t)
    elif ib in (3, -3):
        t = m[0] / math.sqrt(3.0)
        r[:] = t
        if ib > 0:
            r[0, 2] = r[1, 0] = r[2, 1] = -t
        else:
            r[0, 0] = r[1, 1] = r[2, 2] = -t
    elif ib == 4:
        r[0] = (m[0], 0, 0)
        r[1] = (-m[0] / 2.0, m[0] * math.sqrt(3.0) / 2.0, 0)
        r[2] = (0, 0, m[2])
    elif ib == 15:
        r[0] = (m[0], 0, 0)
        r[1] = (m[0] / 2.0, m[0] * math.sqrt(3.0) / 2.0, 0)
        r[2] = (0, 0, m[2])
    elif ib == 7:
        r[0] = (m[0] / 2, m[0] / 2, -m[2] / 2)
        r[1] = (m[0] / 2, -m[0] / 2, m[2] / 2)
        r[2] = (-m[0] / 2, m[0] / 2, m[2] / 2)
    return r


_UNIMODULAR_CACHE: Optional[np.ndarray] = None


def _small_unimodular() -> np.ndarray:
    global _UNIMODULAR_CACHE
    if _UNIMODULAR_CACHE is None:
        allm = np.array(list(itertools.product((-1, 0, 1), repeat=9)), dtype=float)
        allm = allm.reshape(-1, 3, 3)
        det = np.round(np.linalg.det(allm))
        _UNIMODULAR_CACHE = allm[np.abs(det) == 1]
    return _UNIMODULAR_CACHE


@dataclass
class Orientation:
    lattice: np.ndarray        # what to write (= what RMG will use)
    frac: np.ndarray           # fractional coordinates in that lattice
    basis_change: np.ndarray   # integer M with  lattice = M @ input_lattice @ Q
    rotation: np.ndarray       # orthogonal Q
    changed: bool
    mirrored: bool


def orient_for_rmg(lattice: np.ndarray, frac: np.ndarray) -> Orientation:
    """Rewrite a cell so that RMG builds exactly the requested crystal.

    RMG re-orients recognised Bravais lattices to its own convention (``rmg_rotate``). That
    keeps the lattice-vector lengths but may change the angles between the vectors (RMG's
    body-centred cubic vectors differ from QE's ibrav=3) or the orientation (an orthorhombic
    or hexagonal cell given in another orientation), so neither the input fractional nor the
    input Cartesian coordinates are always right in RMG's frame. This function finds the
    integer basis change M (entries -1..1, |det| = 1) and the orthogonal Q with
    ``rmg_rotate(A) = M A Q``, then returns RMG's lattice with fractional coordinates
    ``f M^-1``, i.e. the requested crystal rigidly rotated into RMG's frame. The written
    deck is therefore self-consistent: RMG leaves the written vectors unchanged.

    Body-centred tetragonal cells that RMG detects as such are refused: RMG re-orients them
    with the primitive-vector length standing for both a and c, which is a different lattice.
    """
    a = np.asarray(lattice, dtype=float)
    f = np.asarray(frac, dtype=float).reshape(-1, 3)
    scale = float(np.max(np.abs(a)))
    atol = 1e-10 * max(scale, 1.0)
    if np.allclose(rmg_rotate(a), a, atol=atol):
        return Orientation(a.copy(), f.copy(), np.eye(3), np.eye(3), False, False)
    # RMG's rotation is not always idempotent (its ibrav=3 vectors are re-detected as -3), so
    # the cell written to the deck must be a fixed point of rmg_rotate: RMG then keeps it.
    target = rmg_rotate(a)
    for _ in range(8):
        nxt = rmg_rotate(target)
        if np.allclose(nxt, target, atol=atol):
            break
        target = nxt
    else:  # pragma: no cover - would need a cycle in RMG's rules
        raise ConversionError("RMG's lattice re-orientation does not settle for this cell")
    g_t = target @ target.T
    cands = _small_unimodular()
    ma = cands @ a
    g = ma @ np.transpose(ma, (0, 2, 1))
    tol = 1e-8 * max(scale, 1.0) ** 2
    ok = np.where(np.max(np.abs(g - g_t), axis=(1, 2)) < tol)[0]
    if ok.size == 0:
        raise ConversionError(
            "RMG would re-orient this cell (detected Bravais type %d) to vectors that do not "
            "describe the same lattice; give the cell in RMG's own orientation or use a "
            "conventional cell" % rmg_lat2ibrav(a))
    best = None
    for idx in ok:
        m = cands[idx]
        q = np.linalg.solve(m @ a, target)
        mirrored = np.linalg.det(q) < 0
        cost = (int(mirrored), float(np.abs(m - np.eye(3)).sum()))
        if best is None or cost < best[0]:
            best = (cost, m, q, mirrored)
    _, m, q, mirrored = best
    if not np.allclose(q @ q.T, np.eye(3), atol=1e-8):  # pragma: no cover - guarded by Gram
        raise ConversionError("internal error: lattice re-orientation is not orthogonal")
    f_new = f @ np.linalg.inv(m)
    return Orientation(target, f_new, m, q, True, bool(mirrored))


def apply_orientation(st: Structure, notes: Optional[Notes] = None) -> Tuple[Structure, Orientation]:
    """Return ``st`` expressed in RMG's frame (see :func:`orient_for_rmg`)."""
    o = orient_for_rmg(st.lattice, st.frac)
    if o.changed and notes is not None:
        msg = ("RMG re-orients this Bravais lattice (type %d); the deck holds RMG's lattice "
               "vectors and fractional coordinates for the same crystal, rigidly rotated"
               % rmg_lat2ibrav(st.lattice))
        if not np.allclose(o.basis_change, np.eye(3)):
            msg += ", with a change of primitive basis M = %s" % o.basis_change.astype(int).tolist()
        if o.mirrored:
            msg += " and mirrored (same energies; chirality inverted)"
        notes.info(msg)
    return Structure(o.lattice, list(st.symbols), o.frac, st.movable.copy(), list(st.spin)), o


# ---------------------------------------------------------------------------------------------
# Pseudopotentials
# ---------------------------------------------------------------------------------------------


@dataclass
class UpfHeader:
    path: str
    version: int
    element: Optional[str]
    z_valence: Optional[float]
    pseudo_type: Optional[str]
    has_so: Optional[bool]
    functional: Optional[str]


def read_upf_header(path: str) -> UpfHeader:
    """Read element, valence charge and type from a UPF file (v1 or v2 header)."""
    with open(path, "r", errors="replace") as f:
        text = f.read(200000)
    m = re.search(r"<PP_HEADER\b([^>]*)/?>", text, re.S)
    if m and "=" in m.group(1):
        attrs = dict((k.lower(), v.strip()) for k, v in
                     re.findall(r'(\w+)\s*=\s*"([^"]*)"', m.group(1)))
        z = attrs.get("z_valence")
        so = attrs.get("has_so")
        return UpfHeader(
            path, 2, attrs.get("element", "").strip() or None,
            float(z.replace("D", "E").replace("d", "e")) if z else None,
            attrs.get("pseudo_type"),
            (so.strip().upper().startswith("T") if so else None),
            attrs.get("functional"))
    m = re.search(r"<PP_HEADER>(.*?)</PP_HEADER>", text, re.S)
    if not m:
        raise ValueError("no PP_HEADER found in %s" % path)
    element = z = ptype = func = None
    for line in m.group(1).splitlines():
        low = line.lower()
        tok = line.split()
        if not tok:
            continue
        if "element" in low:
            element = tok[0]
        elif "z valence" in low:
            z = float(tok[0].replace("D", "E").replace("d", "e"))
        elif ptype is None and tok[0].upper() in ("NC", "US", "PAW", "SL"):
            ptype = tok[0].upper()
        elif "exchange-correlation" in low:
            func = " ".join(tok[:4])
    return UpfHeader(path, 1, element, z, ptype, None, func)


def check_upf_for_rmg(h: UpfHeader, notes: Notes) -> None:
    """Report UPF files RMG cannot read (RMG reads UPF v2 XML, NC or US only)."""
    if h.version == 1:
        notes.unsupported(
            "%s is a UPF v1 file; RMG reads UPF v2 (XML) only. Convert it with QE's upfconv.x "
            "-u" % os.path.basename(h.path))
    if h.pseudo_type and h.pseudo_type.upper().startswith("PAW"):
        notes.unsupported("%s is a PAW dataset; RMG supports norm-conserving and ultrasoft "
                          "pseudopotentials only" % os.path.basename(h.path))


_ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As "
    "Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu "
    "Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np "
    "Pu Am Cm Bk Cf Es Fm Md No Lr").split()
ELEMENTS = frozenset(_ELEMENTS)


def element_from_label(label: str) -> Optional[str]:
    """Chemical symbol from a species label such as ``Fe1``, ``Fe_up``, ``Fe_pv``, ``C2``."""
    m = re.match(r"[A-Za-z]+", label.strip())
    if not m:
        return None
    s = m.group(0)
    if len(s) >= 2 and (s[1].islower() or s[:2].isupper()):
        two = s[:2].capitalize()
        if two in ELEMENTS:
            return two
    one = s[0].upper()
    return one if one in ELEMENTS else None


# ---------------------------------------------------------------------------------------------
# k-points
# ---------------------------------------------------------------------------------------------


def check_exx_qmesh(kmesh: Sequence[int], qmesh: Sequence[int]) -> None:
    for k, q in zip(kmesh, qmesh):
        if q < 1 or k % q:
            raise ConversionError(
                "exact-exchange q mesh %s must divide the k mesh %s entry by entry (RMG "
                "XC/Exxbase.cpp)" % (list(qmesh), list(kmesh)))


def write_kpoint_list(deck: RmgDeck, kfrac: np.ndarray, weights: Sequence[float]) -> None:
    """Explicit k-point list in reciprocal-lattice (crystal) units."""
    deck.set("kpoint_mesh", "-1 1 1")
    deck.set("kpoint_units", "Reciprocal lattice")
    deck.set_block("kpoints", ["  %.12f %.12f %.12f %.10g" % (k[0], k[1], k[2], w)
                               for k, w in zip(kfrac, weights)])


# ---------------------------------------------------------------------------------------------
# Small shared utilities
# ---------------------------------------------------------------------------------------------


def nocc_rmg(nel: float, noncollinear: bool = False) -> int:
    """Number of occupied states RMG reserves per k-point (``Misc/init_states.cpp``)."""
    return int(math.ceil(0.5 * nel * (2 if noncollinear else 1) - 1e-9))


def set_scf_thresholds(deck: RmgDeck, etol_ha: float, notes: Notes, hybrid: bool = False,
                       source: str = "conv_thr", spin: bool = False) -> None:
    """Write RMG's SCF (and EXX) stopping criteria for a requested total-energy accuracy etol_ha.

    RMG stops an SCF loop as soon as EITHER ct.scf_accuracy < energy_convergence_criterion OR the
    RMS potential change < rms_convergence_criterion (``RMG/Common/Quench.cpp``, ``Scf.cpp``).
    ct.scf_accuracy is 1/2 int (V_H[rho_out] - V_H[rho_in]) (rho_out - rho_in), the Hartree energy of
    the density residual, which is the quantity QE's conv_thr bounds ("estimated scf accuracy"), so
    without spin the energy criterion is etol_ha itself and the potential criterion,
    1e-3 * sqrt(etol in Ry), is only a backstop.

    The Hartree metric sees the total-charge residual only, not the magnetization residual. With
    spin it can drop below any threshold while the moment is still relaxing (a 2-ZGNR PBE run
    stopped at 4e-13 Ha with RMS[dV] 2e-4 and a moment of 1.00016 instead of 1), so for spin-
    polarized and noncollinear decks the energy criterion is set to RMG's floor (1e-20) and the
    potential criterion carries the convergence. For a hybrid functional the EXX loop stops on
    |exx_delta| < etol_ha (QE: dexx < conv_thr).
    """
    rthr = min(max(1.0e-3 * math.sqrt(2.0 * etol_ha), 1.0e-12), 1.0e-4)
    if spin:
        ethr = 1.0e-20
        notes.info("%s -> rms_convergence_criterion %.3g; energy_convergence_criterion set to RMG's floor "
                   "%.3g because RMG's energy test (Hartree energy of the density residual) does not see "
                   "the magnetization residual" % (source, rthr, ethr))
    else:
        ethr = min(max(etol_ha, 1.0e-20), 1.0e-7)
        notes.info("%s -> energy_convergence_criterion %.3g Ha (the same Hartree-residual measure QE uses) "
                   "and rms_convergence_criterion %.3g: RMG stops on whichever is met first"
                   % (source, ethr, rthr))
    deck.set("energy_convergence_criterion", ethr)
    deck.set("rms_convergence_criterion", rthr)
    if hybrid:
        exx = min(max(etol_ha, 1.0e-14), 1.0e-2)
        deck.set("exx_energy_convergence_criterion", exx)
        notes.info("exx_energy_convergence_criterion %.3g Ha from %s (stop when |exx_delta| is below it); "
                   "older RMG builds ignore it" % (exx, source))


def write_outputs(prog: str, deck: RmgDeck, notes: Notes, out: Optional[str],
                  notes_path: Optional[str], header: str) -> None:
    """Write the deck (file or stdout), the notes sidecar and the notes on stderr."""
    text = deck.render()
    if out and out != "-":
        if os.path.dirname(out):
            os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as f:
            f.write(text)
        npath = notes_path or (out + ".notes")
    else:
        sys.stdout.write(text)
        npath = notes_path
    if npath:
        with open(npath, "w") as f:
            f.write(notes.text(header))
    notes.emit(prog)
