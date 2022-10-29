#!/usr/bin/env python
# calcDelta.py determines the Delta values between two codes. 
# Current version: 3.1
#
# Copyright (C) 2012 Kurt Lejaeghere <Kurt.Lejaeghere@UGent.be>, Center for
# Molecular Modeling (CMM), Ghent University, Ghent, Belgium
#
# calcDelta.py is free software; you can redistribute it and/or modify it 
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your 
# option) any later version.
#
# In addition to the regulations of the GNU Lesser General Public License,
# publications and communications based in parts on this program or on
# parts of this program are required to cite the following articles:
#
# K. Lejaeghere, V. Van Speybroeck, G. Van Oost, and S. Cottenier, "Error 
# estimates for solid-state density-functional theory predictions: an overview
# by means of the ground-state elemental crystals", Crit. Rev. Solid State 39,
# 1-24 (2014).
#
# calcDelta.py is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with calcDelta.py; if not, see <http://www.gnu.org/licenses/>.

# Python and numpy are required to use this script.

from sys import argv, stdout, exit
import os
import getopt
import numpy as np
import sys

if(sys.version[0]==3):
    print('This version of calcDelta requires python version 3')
    exit
    
periodtable = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
               'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
               'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
               'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
               'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
               'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
               'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
               'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
               'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
               'Ba': 56, #'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
#               'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
#               'Er': 68, 'Tm': 69, 'Yb': 70, 
               'Lu': 71, 'Hf': 72, 'Ta': 73,
               'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
               'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, #'At': 85,
               'Rn': 86, #'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
#               'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
#               'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
#               'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
#               'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Uut': 113,
#               'Fl': 114, 'Uup': 115, 'Lv': 116, 'Uus': 117, 'Uuo': 118
               }

# get name of elements, sorted according to atom number
elementlist = [name for (name, num) in sorted(periodtable.items(), 
               key=lambda x:x[1])]

def main():
    usage = '''
calcDelta.py -- Support script to calculate the Delta factor of a code (v3.1)

Use: python calcDelta.py infile [reffile] [-s|--stdout] [-a|--asymm]
     where reffile and infile refer to files containing the element name, V0, 
     B0, and B1 information (V0 in A^3/atom, B0 in GPa, and B1 dimensionless) 
     in columns
     This command calculates the Delta values of the code in infile compared
     to the one in reffile. When reffile is not explicitly given, WIEN2k.txt 
     is used by default.
     Additional output is printed in Delta-out.txt. The option --stdout can be
     used to explicitly print all elements to standard output (on screen)
     instead.
     The current version of the calcDelta script uses a symmetric integration
     between the two codes. For the (deprecated) reference-based integration
     from the Crit. Rev. Solid State article, use the option --asymm.
     Attention: the presence of WIEN2k.txt in the same folder is required for
     full functionality!

python calcDelta.py --help displays the current instructions
    '''

    try:
        opts, args = getopt.gnu_getopt(argv[1:], "hsa", ["help", "stdout", 
            "asymm"])
    except getopt.GetoptError as err:
        print(err)
        print(usage)
        exit(2)

    usestdout = False
    useasymm = False
    for o, a in opts:
        if o in ("-h", "--help"):
            print(usage)
            exit()
        elif o in ("-s", "--stdout"):
            usestdout = True
        elif o in ("-a", "--asymm"):
            useasymm = True
        else:
            assert False, "unhandled option"

    if len(args) < 1:
        print("Error: Please specify a file to read")
        print(usage)
        exit()

    reffile = 'WIEN2k.txt'
    if len(args) == 2:
        reffile = args[1]

    if not os.path.isfile(reffile):
        print("Error: Cannot find reffile %s" % reffile)
        exit(3)

    try:
        data_f = np.loadtxt(args[0], 
            dtype={'names': ('element', 'V0', 'B0', 'BP'),
            'formats': ('U2', float, float, float)})
        data_w = np.loadtxt(reffile, 
            dtype={'names': ('element', 'V0', 'B0', 'BP'),
            'formats': ('U2', float, float, float)})
    except IOError as err:
        print("Error: Cannot read the input files %s and %s: %s" 
            % (args[0], reffile, err))
        exit(4)

    try:
        len(data_f['element'])
    except TypeError:
        print('Error: ' + argv[1] + ': at least two elements required')
        exit()

    eloverlap = list(set(data_f['element']) & set(data_w['element']))

    if not eloverlap:
        print("Error: the input files have no overlapping elements")
        exit()

    Delta, Deltarel, Delta1 = calcDelta(data_f, data_w, eloverlap, useasymm)

    if usestdout:
        outfile = stdout
    else:
        outfile = open('Delta-out.txt', 'w')

    show_results(outfile, args[0], reffile, Delta, Deltarel, Delta1, eloverlap,
                 useasymm)

    outfile.close()


def calcDelta(data_f, data_w, eloverlap, useasymm):
    """
    Calculate the Delta using the data in data_f, data_w on
    element in eloverlap
    """
    v0w = np.zeros(len(eloverlap))
    b0w = np.zeros(len(eloverlap))
    b1w = np.zeros(len(eloverlap))

    v0f = np.zeros(len(eloverlap))
    b0f = np.zeros(len(eloverlap))
    b1f = np.zeros(len(eloverlap))

    elw = list(data_w['element'])
    elf = list(data_f['element'])

    for i in range(len(eloverlap)):
        searchnr = elw.index(eloverlap[i])
        v0w[i] = data_w['V0'][searchnr]
        b0w[i] = data_w['B0'][searchnr] * 10.**9. / 1.602176565e-19 / 10.**30.
        b1w[i] = data_w['BP'][searchnr]

        searchnr = elf.index(eloverlap[i])
        v0f[i] = data_f['V0'][searchnr]
        b0f[i] = data_f['B0'][searchnr] * 10.**9. / 1.602176565e-19 / 10.**30.
        b1f[i] = data_f['BP'][searchnr]

    vref = 30.
    bref = 100. * 10.**9. / 1.602176565e-19 / 10.**30.

    if useasymm:
        Vi = 0.94 * v0w
        Vf = 1.06 * v0w
    else:
        Vi = 0.94 * (v0w + v0f) / 2.
        Vf = 1.06 * (v0w + v0f) / 2.

    a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
    a2f = 9. * v0f**(7./3.) * b0f / 16. * (14. - 3. * b1f)
    a1f = 9. * v0f**(5./3.) * b0f / 16. * (3. * b1f - 16.)
    a0f = 9. * v0f * b0f / 16. * (6. - b1f)

    a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
    a2w = 9. * v0w**(7./3.) * b0w / 16. * (14. - 3. * b1w)
    a1w = 9. * v0w**(5./3.) * b0w / 16. * (3. * b1w - 16.)
    a0w = 9. * v0w * b0w / 16. * (6. - b1w)

    x = [0, 0, 0, 0, 0, 0, 0]

    x[0] = (a0f - a0w)**2
    x[1] = 6. * (a1f - a1w) * (a0f - a0w)
    x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
    x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
    x[4] = -3./5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
    x[5] = -6./7. * (a3f - a3w) * (a2f - a2w)
    x[6] = -1./3. * (a3f - a3w)**2.

    y = [0, 0, 0, 0, 0, 0, 0]

    y[0] = (a0f + a0w)**2 / 4.
    y[1] = 3. * (a1f + a1w) * (a0f + a0w) / 2.
    y[2] = -3. * (2. * (a2f + a2w) * (a0f + a0w) + (a1f + a1w)**2.) / 4.
    y[3] = -(a3f + a3w) * (a0f + a0w) / 2. - (a2f + a2w) * (a1f + a1w) / 2.
    y[4] = -3./20. * (2. * (a3f + a3w) * (a1f + a1w) + (a2f + a2w)**2.)
    y[5] = -3./14. * (a3f + a3w) * (a2f + a2w)
    y[6] = -1./12. * (a3f + a3w)**2.

    Fi = np.zeros_like(Vi)
    Ff = np.zeros_like(Vf)

    Gi = np.zeros_like(Vi)
    Gf = np.zeros_like(Vf)

    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2.*n-3.)/3.)
        Ff = Ff + x[n] * Vf**(-(2.*n-3.)/3.)

        Gi = Gi + y[n] * Vi**(-(2.*n-3.)/3.)
        Gf = Gf + y[n] * Vf**(-(2.*n-3.)/3.)

    Delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))
    Deltarel = 100. * np.sqrt((Ff - Fi) / (Gf - Gi))
    if useasymm:
        Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
                 / v0w / b0w * vref * bref
    else: 
        Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
                 / (v0w + v0f) / (b0w + b0f) * 4. * vref * bref

    return Delta, Deltarel, Delta1


def show_results(outfile, infile, reffile, Delta, Deltarel, Delta1, eloverlap,
                 useasymm):
    """
    Print the result to a file descriptor
    """
    Dmax = Delta.argmax()
    Dmin = Delta.argmin()
    Dmean = Delta.mean()
    Dstdev = Delta.std()

    Drelmax = Deltarel.argmax()
    Drelmin = Deltarel.argmin()
    Deltarelav = Deltarel.mean()
    Drelstdev = Deltarel.std()

    D1max = Delta1.argmax()
    D1min = Delta1.argmin()
    Delta1av = Delta1.mean()
    D1stdev = Delta1.std()

    total = len(eloverlap)

    outfile.write('--------------------\n')
    outfile.write('# Delta values of ' + infile + ' with respect to ' + 
        reffile + ' (in meV/atom)\n')
    outfile.write('# (%i elements of %i included)\n' 
        % (total, len(elementlist)))
    outfile.write('# calculated with calcDelta.py version 3.1 ')
    if useasymm:
        outfile.write('(asymmetric mode) \n')
    else:
        outfile.write('\n')
    outfile.write('# from left to right: Delta [meV/atom] - relative Delta [%]'
        + ' - Delta1 [meV/atom]\n')
    outfile.write('--------------------\n')

    for el in elementlist:
        while True:
            try:
                i = eloverlap.index(el)
                outfile.write(eloverlap[i] + '\t %.3f\t %.1f\t%.3f \n' 
                    % (Delta[i], Deltarel[i], Delta1[i]))
                break
            except ValueError:
                outfile.write(el + '\t N/A \t N/A \tN/A \n')
            break

    outfile.write('--------------------\n')
    outfile.write('np.mean  %.3f\t %.1f\t%.3f\n' 
        % (Dmean, Deltarelav, Delta1av))
    outfile.write('np.std   %.3f\t %.1f\t%.3f\n' 
        % (Dstdev, Drelstdev, D1stdev))
    outfile.write('np.max   %.3f\t %.1f\t%.3f \t (%.2s, %.2s, %.2s)\n'
        % (Delta[Dmax], Deltarel[Drelmax], Delta1[D1max], eloverlap[Dmax],
        eloverlap[Drelmax], eloverlap[D1max]))
    outfile.write('np.min   %.3f\t %.1f\t%.3f \t (%.2s, %.2s, %.2s)\n'
        % (Delta[Dmin], Deltarel[Drelmin], Delta1[D1min], eloverlap[Dmin],
        eloverlap[Drelmin], eloverlap[D1min]))
    outfile.write('--------------------\n')

if __name__ == "__main__":
    main()
