#!/usr/bin/python3

import numpy
import os
import stat
import math
import pymatgen
from pymatgen.io.cif import CifParser
from pymatgen.core import lattice
from pymatgen.core import structure
from pymatgen.io.pwscf import PWInput, PWInputError, PWOutput
import re

# Volumes ratios for delta test
volume_lists = [0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03]
#pseudo_extra = ["Ar", "Bi", "He", "Kr", "Lu", "Ne", "Po", "Rn", "Xe"]
pseudo_extra = []


semicore_list =["Sn", "Tl"]
FM_list = ["Fe", "Co", "Ni"]
AFM_list = ["Cr", "Mn", "O"]

grid_spacing = 0.12  # in unit of Angstrom
k_delta = 0.2  # in unit of (Anstrom^-1)
k_parall = 3;
#pp = 'davidson_broad0.01eV_evengrid_MP'
pp = 'test'


def main():

    jobline = """
#!/bin/bash
cp ~/bin/rmg-cpu .
"""
    line_run = 'mpirun -np 24 ../../rmg-cpu input\n'
    shell_script_line=''

    dirname = 'Delta_benchmark_dk_%4.2f_hx_%4.2f_%s'%(k_delta, grid_spacing,pp)
    if not os.path.exists(dirname): os.mkdir(dirname)
    os.chdir(dirname)

    for filename in os.listdir('../../CIFs/'):
        species = filename.split('.')[0]
        parser = CifParser('../../CIFs/' + filename)
        struc_conv = parser.get_structures(0)[0]
        struc_prim = parser.get_structures(1)[0]
        structure = struc_prim
        
        # RMG has fast operators for orthogonal cells, hex, fcc and bcc so we try to use
        # the smallest primitive cell that matches one of those if possible
        a_p = struc_prim.lattice.a
        b_p = struc_prim.lattice.b
        c_p = struc_prim.lattice.c
        a_c = struc_conv.lattice.a
        b_c = struc_conv.lattice.b
        c_c = struc_conv.lattice.c
        alpha_p = struc_prim.lattice.alpha
        beta_p = struc_prim.lattice.beta
        gamma_p = struc_prim.lattice.gamma
        alpha_c = struc_conv.lattice.alpha
        beta_c = struc_conv.lattice.beta
        gamma_c = struc_conv.lattice.gamma

        c_is_ortho = struc_conv.lattice.is_orthogonal
        c_is_hex = struc_conv.lattice.is_hexagonal
        p_is_ortho = struc_prim.lattice.is_orthogonal
        p_is_hex = struc_prim.lattice.is_hexagonal
       
        
        if p_is_ortho:
            structure = struc_prim
        elif p_is_hex:
            structure = struc_prim
        elif c_is_ortho:
            structure = struc_conv
        elif c_is_hex:
            structure = struc_conv



        # Base volume
        v0 = structure.volume
        #print(structure.lattice)
        rmesh = [0, 0, 0]
        vectors = structure.lattice.matrix
        i = 0
        for v in vectors:
            #vmags[i] = numpy.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
            #print(vmags[i])
            i = i + 1

        if not os.path.exists(species): os.mkdir(species)
        os.chdir(species)

        shell_script_line += 'cd %s\n'%(species)
        shell_script_line += """grep "volume and e" */*.log |awk '{print $8 "  " $9}'|sort >vol_energy.dat"""

        shell_script_line += '\n'
        shell_script_line += 'python ../../../eosfit.py vol_energy.dat\n'
        shell_script_line += 'if [ -f "vol_energy.dat.eosout" ]\n'
        shell_script_line += 'then\n'
        shell_script_line += 'echo "%s" >>../tem.txt\n'%(species)
        shell_script_line += 'awk "NR==3" vol_energy.dat.eosout>>../tem.txt\n'
        shell_script_line += 'fi\n'
        shell_script_line += 'cd ..\n'


        input_for_rmg(species, structure, grid_spacing, v0)

        os.chdir('..')
        jobline += 'cd '+species+'/volume_0.97\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_0.98\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_0.99\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.0\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.01\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.02\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.03\n'
        jobline += line_run
        jobline += 'cd ../..\n'

    job_blue_water = open('job_blue_water', 'w')
    job_blue_water.write(jobline)
    job_blue_water.close()

    shell_script_line += """paste - - -d"    " <tem.txt >RMG.txt\n"""
    shell_script_line += 'python ../../calcDelta.py RMG.txt ../../WIEN2k.txt\n'
    f = open('eosfit.sh', 'w')
    f.write('rm tem.txt\n')
    f.write(shell_script_line)
    f.close()
    os.chmod('eosfit.sh', stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR)

def input_for_rmg(species, structure, grid_spacing, v0):
    default_line = """
crds_units = "Angstrom"
lattice_units = "Angstrom"
atomic_coordinate_type = "Cell Relative"
potential_grid_refinement="2"
occupations_type = "MethfesselPaxton"
occupation_electron_temperature_eV = "1.0e-2"
charge_density_mixing = "0.15"
potential_acceleration_constant_step = "0.0"
charge_mixing_type = "Broyden"
output_wave_function_file = "/dev/null"
ecutwfc="140.0"
rms_convergence_criterion = "1e-7"
start_mode="LCAO Start"
max_scf_steps = "100"
calculation_mode="Quench Electrons"
localize_projectors = "false"
localize_localpp = "false"
kohn_sham_solver="davidson"
kpoint_distribution = "%d"
"""%(k_parall)

    semi_core_line = """
kohn_sham_mg_timestep = "0.1"
kohn_sham_time_step = "0.3"
kohn_sham_coarse_time_step = "0.1"
"""

    atom_format = " %s     %.12e    %.12e    %.12e      1  %f\n"

    atom_list = []
    for site in structure:
        x = site.frac_coords[0]
        y = site.frac_coords[1]
        z = site.frac_coords[2]
        atom_list.append([species, x, y, z])

    final_atom_list = []
    final_atom_list = atom_list
    _positions_line = ''
    mag = 0.0;
    sign = 1
    if(species in FM_list):
        mag = 0.25
        sign = 1
        _positions_line +='spin_polarization="true"\n'

    if(species in AFM_list):
        mag = 0.25
        sign = -1
        _positions_line +='spin_polarization="true"\n'
    _positions_line += 'atoms=\n"\n'
    for atom in final_atom_list:
        _positions_line += ( atom_format % (atom[0], atom[1], atom[2], atom[3], mag) )
        mag = mag * sign

    _positions_line += '"\n'

    a0 = structure.lattice.a
    b0 = structure.lattice.b
    c0 = structure.lattice.c

    (nx, ny, nz) = SetupGrid(a0, b0, c0, grid_spacing)
    kx = int(2.0 * 3.1415926/a0/k_delta)
    ky = int(2.0 * 3.1415926/b0/k_delta)
    kz = int(2.0 * 3.1415926/c0/k_delta)
    if (kx == 0): kx = 1
    if (ky == 0): ky = 1
    if (kz == 0): kz = 1
    _positions_line += 'kpoint_mesh = "%d %d %d"\n'%(kx, ky, kz)
    _positions_line += 'wavefunction_grid = "%d %d %d"\n'%(nx, ny, nz)

    jobfile = open('job.bat', 'w')
    jobfile.write('rm */*.log */*.options\n')
    jobrun = 'mpirun -n 16 ~/SVN/codes/rmg-cpu input\n'

    pp_line='\n'
    if(pp == 'ONCV'):
        pp_line = ' pseudopotential = "%s  ../../../ONCV/%s_ONCV_PBE-1.0.upf"\n'%(atom_list[0][0], atom_list[0][0])

    if(species in pseudo_extra):
        pp_line ='pseudopotential = "%s ../../../PBE_MT_FHI/%s.pbe-mt_fhi.txt"\n '%(species, species)

    semi_core_line ='\n'
    if (species in semicore_list):
        semi_core_line = """
kohn_sham_mg_timestep = "0.1"
kohn_sham_time_step = "0.3"
kohn_sham_coarse_time_step = "0.1"
"""

    for vol in volume_lists:
        dir_name = 'volume_' + str(vol)
        jobfile.write('cd ' + dir_name + '\n')
        jobfile.write(jobrun)
        jobfile.write('cd ..\n')
        if not os.path.exists(dir_name): os.mkdir(dir_name)
        f = open(dir_name+'/input', 'w')
        structure.scale_lattice(vol * v0)
        vectors = structure.lattice.matrix
        lattice_line = "lattice_vector = \""
        for v in vectors:
            lattice_line += "\n%14.8f    %14.8f    %14.8f" % (v[0],v[1],v[2])
        lattice_line += "\n\"\n\n"

        f.write(lattice_line)
        f.write(_positions_line)
        f.write(default_line);
        f.write(pp_line)
        f.write(semi_core_line)
        f.close()

    jobfile.close()
    os.chmod('job.bat', stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR)



def SetupGrid(a0, b0, c0, grid_spacing):
    nx = int(a0/grid_spacing)
    ny = int(b0/grid_spacing)
    nz = int(c0/grid_spacing)
#    nx = int ((nx + 7)/8) * 8
#    ny = int((ny + 7)/8) * 8
#    nz = int((nz + 7)/8) * 8

    nx = int ((nx + 1)/2) * 2
    ny = int((ny + 1)/2) * 2
    nz = int((nz + 1)/2) * 2


    spacing_list =[a0/nx, b0/ny, c0/nz]
    max_sp = max(spacing_list)
    min_sp = min(spacing_list)
    anisotropy = max_sp/min_sp
    max_iter = 0
    while(anisotropy >1.1 and max_iter < 10):
        max_iter +=1
        nx += 2
        xspacing = a0/nx
        ny = b0/xspacing
        ny = int((ny + 1)/2) * 2
        nz = c0/xspacing
        nz = int((nz + 1)/2) * 2
        spacing_list =[a0/nx, b0/ny, c0/nz]
        max_sp = max(spacing_list)
        min_sp = min(spacing_list)
        anisotropy = max_sp/min_sp
        print (anisotropy)


    return(nx, ny, nz)


if __name__ == '__main__':
    main()

