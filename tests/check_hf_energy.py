#!/usr/bin/env python3
"""Compare the post-SCF Hartree-Fock energy printed for QMCPACK with hf_reference.txt.

hf_reference.txt: "<energy in Ha> <tolerance in Ha>". Prints "pass" or "fail" for CTest.
"""
import glob
import re
import sys

ref, tol = (float(x) for x in open('hf_reference.txt').read().split()[:2])
logs = sorted(glob.glob('input.*.log'))
if not logs:
    print('fail: no log file')
    sys.exit(1)
values = re.findall(r'total energy Hartree Fock\(\+Madelung\)\s*=\s*(-?[0-9.]+)', open(logs[-1]).read())
if not values:
    print('fail: no Hartree Fock energy in %s' % logs[-1])
    sys.exit(1)
e = float(values[-1])
status = 'pass' if abs(e - ref) <= tol else 'fail'
print('Hartree Fock(+Madelung) energy %.8f Ha, reference %.8f Ha, tolerance %.1e: %s' % (e, ref, tol, status))
for line in open('test_energy.log'):
    print(line, end='')
