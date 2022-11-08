import os
import sys
import string
import copy
import math
from datetime import datetime
from optparse import OptionParser, OptionGroup
import warnings
import CifFile
import subprocess
from utils import *
from uctools import *

def token_to_value(all_lines, token):
    for line in all_lines:
        if(token in line):
            t_pos = line.find(token)
            v_pos_b = line.find("=", t_pos)
            v_pos_e = line.find(",", t_pos)
            return line[v_pos_b+1:v_pos_e]
    return None    
class rmg_interface():

    def cif2cell(self, cif_file=None): 
        #################################################################
        # Open and read CIF file
        cf = CifFile.ReadCif(cif_file)

        self.atom_unit = "Absolute"
        ##############################################
        # Get blocks
        cfkeys = list(cf.keys())
        cb = cf.get(cfkeys[0])
        # Get reference data
        ref = ReferenceData()
        ref.getFromCIF(cb)
        # Get cell data
        cd = CellData()
        cd.force = True
        cd.getFromCIF(cb)


        ##############################################
        # Generate cell
        if self.reducetoprim:
            cd.primitive()
        else:
            cd.conventional()

        inputcell = copy.copy(cd)

        self.cell = cd

        self.cell.newunit("angstrom")

        self.ibrav = 0
        t = LatticeMatrix(self.cell.latticevectors)
        abc_length = [0,0,0]
        for i in range(3):
            for j in range(3):
                t[i][j] = self.cell.latticevectors[i][j]*self.cell.lengthscale
                abc_length[i] += t[i][j] * t[i][j]
            abc_length[i] = math.sqrt(abc_length[i])    
        self.cell.a = abc_length[0]
        self.cell.b = abc_length[1]
        self.cell.c = abc_length[2]

        ortho = abs(self.cell.a - t[0][0]) < 1.0e-5 
        ortho = ortho and abs(self.cell.b - t[1][1]) < 1.0e-5
        ortho = ortho and abs(self.cell.c - t[2][2]) < 1.0e-5

        if ortho: self.ibrav = 8
        system = self.cell.crystal_system()
        setting = self.cell.spacegroupsetting
        if system == 'cubic':
            if self.cell.primcell:
                if setting == 'P':
                    self.ibrav = 1
                elif setting == 'F':
                    self.ibrav = 2
                elif setting == 'I':
                    self.ibrav = 3
            else:
                self.ibrav = 1
        if system == 'hexagonal':
            if self.cell.primcell:
                if setting == 'P':
                    self.ibrav = 4
            #elif setting == 'R':
            #    return 5
        if system == 'tetragonal':
            if self.cell.primcell:
                if setting == 'P':
                    self.ibrav = 8
                elif setting == 'I':
                    self.ibrav = 7
            else:
                self.ibrav = 8

        self.atoms = []
        for a in self.cell.atomdata:
            for b in a:
                self.atoms.append([b.spcstring(), b.position[0],b.position[1],b.position[2]])
        #print(self.atoms)

    def cell2rmg(self, mag):
        filestring = "#****  LATTICE and ATOMS  ****   \n"
        #
        # some default input options
        brav_type = {
            0:"None",
            1:"Cubic Primitive",
            2:"Cubic Face Centered",
            3:"Cubic Body Centered",
            4:"Hexagonal Primitive",
            8:"Orthorhombic Primitive"
        }
        filestring += 'bravais_lattice_type="%s"  \n'%brav_type[self.ibrav]
        if self.cell.unit == "angstrom":
            filestring += 'crds_units = "Angstrom"  \n'
            filestring += 'lattice_units = "Angstrom"  \n'
        elif self.cell.unit == "bohr":
            filestring += 'crds_units = "Bohr"  \n'
            filestring += 'lattice_units = "Bohr"  \n'
        else:
            st.markdown("WARNING: unit = ",self.cell.unit)

        t = LatticeMatrix(self.cell.latticevectors)
        if self.ibrav !=0:
            filestring += 'a_length="%16.8f"  \n'%self.cell.a
            filestring += 'b_length="%16.8f"  \n'%self.cell.b
            filestring += 'c_length="%16.8f"  \n'%self.cell.c
        else:
            t = LatticeMatrix(self.cell.latticevectors)
            filestring += 'lattice_vector="  \n'
            for i in range(3):
                for j in range(3):
                    filestring += " %.12e "%(self.cell.latticevectors[i][j] * self.cell.lengthscale)
                filestring += '  \n'    
            filestring += '"  \n'

        filestring += 'atomic_coordinate_type = "%s"  \n'%self.atom_unit 
        filestring += 'atoms="  \n'
        atom_format = "%s  %.12e %.12e %.12e"
        iatom = 0
        for a in self.atoms:
            b = a[0]
            if b[len(b) -1].isdigit():
                b = b[:len(b)-1]

            filestring += atom_format%(b,a[1], a[2], a[3])
            filestring += "  1 1 1 %6.2f %6.2f %6.2f  \n"%(mag[iatom][0],mag[iatom][1],mag[iatom][2])
            iatom += 1
        filestring += '"  \n'

        return filestring

    def __init__(self, filename, filetype):
        self.reducetoprim = True
        if not os.path.exists(filename):
            sys.stderr.write("***Error: The file "+filename+" could not be found.\n")
            sys.exit(2)
        if filetype == "cif":
            self.cif2cell(filename)
        else:
            sys.stderr.write("***Error: The file "+filename+" is not .cif\n")
            sys.exit(2)

        # set up species list
        tmp = set([])
        tmp_AFM = set([])
        for a in self.atoms:
            b = a[0]
            if b[len(b) -1].isdigit():
                b = b[:len(b)-1]
            tmp_AFM.add(a[0])
            tmp.add(b)
        self.species = list(tmp)
        self.species_AFM = list(tmp_AFM)
        #self.rmginput = self.cell2rmg(mag)
