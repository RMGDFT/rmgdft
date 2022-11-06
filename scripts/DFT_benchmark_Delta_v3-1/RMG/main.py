#!/usr/bin/python3
#
# RMG specific setup script for delta test. Job run files must be configured for a specific system.



#from numpy import *
import numpy  
import os
import stat
import math
import tarfile


from numpy import sqrt, arccos, fabs, pi, cos, sin
from rmg_parser import *

pi = 3.141592653589793238

# Adjust grids here
grids = {
"Ag":"""wavefunction_grid="20 20 20 """,
"Al":"""wavefunction_grid="20 20 20 """,
"Ar":"""wavefunction_grid="32 32 32 """,
"As":"""wavefunction_grid="28 28 76 """,
"Au":"""wavefunction_grid="24 24 24 """,
"Ba":"""wavefunction_grid="36 36 36 """,
"Be":"""wavefunction_grid="16 16 24 """,
"Bi":"""wavefunction_grid="32 32 88 """,
"Br":"""wavefunction_grid="48 24 52 """,
"B":"""wavefunction_grid="44 44 44 """,
"Ca":"""wavefunction_grid="28 28 28 """,
"Cd":"""wavefunction_grid="20 20 40 """,
"Cl":"""wavefunction_grid="44 24 52 """,
"Co":"""wavefunction_grid="16 16 28 """,
"Cr":"""wavefunction_grid="16 16 16 """,
"Cs":"""wavefunction_grid="44 44 44 """,
"Cu":"""wavefunction_grid="20 20 20 """,
"C":"""wavefunction_grid="16 16 56 """,
"Fe":"""wavefunction_grid="20 20 20 """,
"F":"""wavefunction_grid="44 48 24 """,
"Ga":"""wavefunction_grid="32 56 32 """,
"Ge":"""wavefunction_grid="28 28 28 """,
"He":"""wavefunction_grid="24 24 44 """,
"Hf":"""wavefunction_grid="24 24 44 """,
"Hg":"""wavefunction_grid="24 24 20 """,
"H":"""wavefunction_grid="28 28 36 """,
"In":"""wavefunction_grid="24 24 36 """,
"Ir":"""wavefunction_grid="24 24 24 """,
"I":"""wavefunction_grid="48 24 56 """,
"Kr":"""wavefunction_grid="40 40 40 """,
"K":"""wavefunction_grid="36 36 36 """,
"Li":"""wavefunction_grid="52 52 52 """,
"Lu":"""wavefunction_grid="36 36 64 """,
"Mg":"""wavefunction_grid="24 24 36 """,
"Mn":"""wavefunction_grid="24 24 36 """,
"Mo":"""wavefunction_grid="20 20 20 """,
"Na":"""wavefunction_grid="64 64 64 """,
"Nb":"""wavefunction_grid="24 24 24 """,
"Ne":"""wavefunction_grid="24 24 24 """,
"Ni":"""wavefunction_grid="20 20 20 """,
"N":"""wavefunction_grid="48 48 48 """,
"Os":"""wavefunction_grid="24 24 40 """,
"O":"""wavefunction_grid="32 28 28 """,
"Pb":"""wavefunction_grid="24 24 24 """,
"Pd":"""wavefunction_grid="20 20 20 """,
"Po":"""wavefunction_grid="24 24 24 """,
"Pt":"""wavefunction_grid="24 24 24 """,
"P":"""wavefunction_grid="24 80 32 """,
"Rb":"""wavefunction_grid="40 40 40 """,
"Re":"""wavefunction_grid="16 16 24 """,
"Rh":"""wavefunction_grid="20 20 20 """,
"Rn":"""wavefunction_grid="40 40 40 """,
"Ru":"""wavefunction_grid="24 24 36 """,
"Sb":"""wavefunction_grid="24 24 64 """,
"Sc":"""wavefunction_grid="24 24 36 """,
"Se":"""wavefunction_grid="32 32 36 """,
"Si":"""wavefunction_grid="20 20 20 """,
"Sn":"""wavefunction_grid="32 32 32 """,
"Sr":"""wavefunction_grid="20 20 20 """,
"S":"""wavefunction_grid="20 20 20 """,
"Ta":"""wavefunction_grid="24 24 24 """,
"Tc":"""wavefunction_grid="20 20 32 """,
"Te":"""wavefunction_grid="32 32 40 """,
"Ti":"""wavefunction_grid="20 20 32 """,
"Tl":"""wavefunction_grid="24 24 40 """,
"V":"""wavefunction_grid="20 20 20 """,
"W":"""wavefunction_grid="16 16 16 """,
"Xe":"""wavefunction_grid="36 36 36 """,
"Y":"""wavefunction_grid="24 24 40 """,
"Zn":"""wavefunction_grid="20 20 40""",
"Zr":"""wavefunction_grid="24 24 36 """}

#external pps here
external_pps = {
"Mn":"""Mn.pbe-spn-rrkjus_psl.0.3.1.UPF""",
"Pt":"""Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF""",
"Os":"""Os.pbe-spfn-rrkjus_psl.1.0.0.UPF""",
"Lu":"""Lu.pbe-spdfn-rrkjus_psl.1.0.0.UPF""",
"Ir":"""Ir.pbe-spfn-rrkjus_psl.1.0.0.UPF""",
"Si":"""Si.pbe-n-rrkjus_psl.1.0.0.UPF""",
"Ga":"""Ga.pbe-dn-rrkjus_psl.1.0.0.UPF""",
"Po":"""Po.pbe-dn-rrkjus_psl.1.0.0.UPF""",
"As":"""As.pbe-n-rrkjus_psl.1.0.0.UPF""",
"N":"""N.oncvpsp.upf""",
"He":"""He_ONCV_PBE-1.0.oncvpsp.upf""",
"Na":"""Na_ONCV_PBE-1.0.oncvpsp.upf""",
"Zn":"""Zn.pbe-spn-rrkjus_psl.1.0.0.UPF"""}

#extra options here
extra_opts = {
"Bi":"""potential_grid_refinement="3""",
"Lu":"""potential_grid_refinement="3""",
"Nb":"""potential_grid_refinement="3"\nkpoint_mesh = "23 23 23 """,
"Ta":"""kpoint_mesh = "31 31 31 """,
"Pb":"""kpoint_mesh = "25 25 25 """,
"Ge":"""kpoint_mesh = "13 13 13 """,
"K":"""kpoint_mesh = "21 21 21 """
}

# We use the pslibrary pseudopotentials with some atomic species
# with info given below. The script does not yet set these up automatically
#
# Mn uses Mn.pbe-spn-rrkjus_psl.0.3.1.UPF with a "24 24 36" grid
# Pt uses Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF with a "24 24 24" grid
# Os uses Os.pbe-spfn-rrkjus_psl.1.0.0.UPF with a "24 24 40" grid
# Lu uses Lu.pbe-spdfn-rrkjus_psl.1.0.0.UPF with a "36 36 64" grid since it has f-electrons in valence
# Ir uses Ir.pbe-spfn-rrkjus_psl.1.0.0.UPF with a "24 24 24" grid and pgr=3
# Si uses Si.pbe-n-rrkjus_psl.1.0.0.UPF with a "20 20 20" grid
# Ga uses Ga.pbe-dn-rrkjus_psl.1.0.0.UPF with a "32 56 32" grid
# Po uses Po.pbe-dn-rrkjus_psl.1.0.0.UPF with a ""24 24 24" grid and pgr=3 and fd_order=10
# As uses As.pbe-n-rrkjus_psl.1.0.0.UPF with a "28 28 76" grid
# Zn uses Zn.pbe-spn-rrkjus_psl.1.0.0.UPF with a "20 20 40" grid
# The script does not automaticially set these up yet.
#
# Ag uses uspp but needs wavefunction_grid = "20 20 20"
# Ni uses uspp but needs wavefunction_grid = "20 20 20"
# Rn uses nc_accuracy but needs wavefunction_grid= "40 40 40"
# Fe uses sg15 but needs wavefunction_grid = "20 20 20"
# Au uses sg15 but needs wavefunction_grid = "24 24 24"
# The base value is used for ncpp
grid_spacing_base = 0.14 # in unit of Angstrom
volume_lists = [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06]
#pseudo_extra = ["Ar", "Bi", "He", "Kr", "Lu", "Ne", "Po", "Rn", "Xe"]

# Rmg default pseudo is SG15. If a species is in this list then USPP will be used.
pseudo_uspp=["Ag","Ni","Sb","W","C","Re","Cr","I","Br","Cl","Hg","B","S","Zn","Ir"]

# ONCV with core corrections
#pseudo_nc = ["Ar", "Kr", "Rn", "W","V"];
# Cu needs a 20x20x20 grid
pseudo_nc=["Fe","Nb","V","Bi","Po","Rn","Kr"]

# Use semilocal form (only available with sg15)
pseudo_sl=["Ar"]

# Extra cutoff required
high_cutoff = ["Cu","Ni","Kr","Pt","B","S", "Zn"]
high_high_cutoff = ["Ir"]

# Higher temperature for convergence
high_temp = ["Cr","Na"]

# Default solver is davidson but multigrid will be used if a species is in this list
run_mg = ["O"];

# Temporary for testing use fft fd
#fft_fd = ["B","Li","O","S","Na","F"]

FM_list = ["Fe", "Co", "Ni"]
# spin up, spin down
AFM_list1 = ["Cr", "Mn"]
#spin up, up, down, down
AFM_list2 = ["O"]

# Denser kpoint mesh for elements in this list. Mostly fcc metals
high_k_list = ["Au", "Pt", "Rh", "Ag", "Ir","Cu","Rn","Cs","Ar"]

# even denser kpoint mesh for elements in this list
#veryhigh_k = ["Au", "Pt"];

# list of species use CIFs instead of primCIFs
# I, P, Br, Cl, Ga, In, Hg: use CIFs will have orthorhombic 
# primCIFs will have monoclinic
#Sb,As, Na, S, B, Li, Bi symmetry not recongnized by cif2cell script
use_CIFs = ["I","P","Br","Cl","Ga","In","Hg","Lu","Sb","As","Na", "S","B","Li","Bi","F"]

k_delta = 0.08  # in unit of (Anstrom^-1)
k_parall = 8;
#pp = 'MIXED'
#pp = 'SG15'
pp='test'
#pp = 'ONCV'

def main():
    jobline = """
#!/bin/bash
cp ~/bin/rmg-cpu .
"""

    jobline_desktop = ''
#    line_run = 'aprun -n 128 -N 8 -d 4 -cc numa_node ../../rmg-cpu input&\n'
    #line_run = 'srun -AMAT189_crusher --ntasks=32 -u -c8 --gpus-per-node=8 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=single:1 --cpu-bind=sockets  ../../rmg-cpu-crusher input\n'
    line_run = 'srun -AMAT189_crusher --ntasks=256 -u --gpus-per-node=8 --ntasks-per-node=64 --cpu-bind=sockets  ../../rmg-cpu-crusher input\n'
#    line_run = 'mpirun -np 24 --bind-to core ~/bin/rmg-cpu input\n'
    
    shell_script_line=''
    dirname = 'Delta_benchmark_dk_%4.2f_hx_%4.2f_%s'%(k_delta, grid_spacing_base, pp)
    if not os.path.exists(dirname): os.mkdir(dirname)
    os.chdir(dirname)
    print("dir", dirname)
    for filename in os.listdir('../../primCIFs/'):
        species = filename.split('.')[0]
    
        print("starting spec", species)
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

        grid_spacing = grid_spacing_base;
        if(species in pseudo_uspp):
            grid_spacing = 1.25*grid_spacing_base;

        if(species in high_cutoff):
            grid_spacing = grid_spacing_base/1.25;

        if(species in high_cutoff and species in pseudo_uspp):
            grid_spacing = grid_spacing_base/1.2;

        if(species in high_high_cutoff):
            grid_spacing = grid_spacing_base/1.75;

        if species in use_CIFs:
            input_for_rmg(species, '../../../CIFs/'+filename,  grid_spacing)
        else:
            input_for_rmg(species, '../../../primCIFs/'+filename,  grid_spacing)

        os.chdir('..')
        jobline_desktop += 'cd ' + species +'\n'
        jobline_desktop += './job.bat\n'
        jobline_desktop += 'cd ..\n'

        jobline += 'cd '+species+'/volume_0.94\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_0.96\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_0.98\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.00\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.02\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.04\n'
        jobline += line_run
        jobline += 'cd ../../'+species+'/volume_1.06\n'
        jobline += line_run
        jobline += 'cd ../..\n'


    job_blue_water = open('job_blue_water', 'w')
    job_blue_water.write(jobline)
    job_blue_water.write('wait\n')
    job_blue_water.close()
    
    shell_script_line += """paste - - -d"    " <tem.txt >RMG.txt\n"""
    shell_script_line += 'python ../../calcDelta.py RMG.txt ../../WIEN2k.txt\n'
    f = open('eosfit.sh', 'w')
    f.write('rm tem.txt\n')
    f.write(shell_script_line)
    f.close()
    os.chmod('eosfit.sh', stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR)

    f = open('job_desktop.sh', 'w')
    f.write(jobline_desktop)
    f.close()
    os.chmod('job_desktop.sh', stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR)

    os.chdir('..')
    #os.system("tar cvf delta.tar " + dirname + ' ONCV PBE_MT_FHI')

def input_for_rmg(species, ciffile, grid_spacing):



    crmg = rmg_interface(ciffile, "cif")

    default_line = """
crds_units = "Angstrom"
lattice_units = "Angstrom"
atomic_coordinate_type = "Cell Relative"
occupations_type = "MethfesselPaxton"
output_wave_function_file = "/dev/null"
rms_convergence_criterion = "1e-7"
start_mode="LCAO Start"
max_scf_steps = "100"
calculation_mode="Quench Electrons"
localize_projectors = "false"
localize_localpp = "false"
kpoint_distribution = "%d"
verbose="true"
energy_convergence_criterion = "1.00000000e-9"
"""%(k_parall)

    atom_format = " %s     %.12e    %.12e    %.12e      1  %f\n"

    _positions_line = ''
    mag = 0.0;
    sign = []
    for atom in crmg.atoms:
        sign.append(1);
    if(species in FM_list):
        mag = 0.25
        sign = [1,1]
        _positions_line +='spin_polarization="true"\n'

    if(species in AFM_list1):
        mag = 0.25
        sign = [1,-1]
        _positions_line +='spin_polarization="true"\n'
        _positions_line +='AFM="true"\n'
    if(species in AFM_list2):
        mag = 0.25
        sign = [1,1,-1,-1]
        _positions_line +='spin_polarization="true"\n'
        _positions_line +='AFM="true"\n'
    _positions_line += 'atoms=\n"\n'
    iatom = 0
    for atom in crmg.atoms:
        mag1 = mag * sign[iatom]
        iatom += 1
        _positions_line += ( atom_format % (atom[0], atom[1], atom[2], atom[3], mag1) )

    _positions_line += '"\n'

    a0 = crmg.cell.a
    b0 = crmg.cell.b
    c0 = crmg.cell.c

# Adjust grid spacing based on ibrav type
    print("ibrav", crmg.ibrav)
    if(crmg.ibrav == 2):
        grid_spacing = grid_spacing * sqrt(2.0);
#    if(crmg.ibrav == 3):
#        grid_spacing = grid_spacing*sqrt(3.0/2.0);
    if(crmg.ibrav == 4):
        grid_spacing = grid_spacing;
    
    
    (nx, ny, nz) = SetupGrid(a0, b0, c0, grid_spacing)
    kx = int(2.0 * 3.1415926/a0/k_delta)
    ky = int(2.0 * 3.1415926/b0/k_delta)
    kz = int(2.0 * 3.1415926/c0/k_delta)
    if (kx == 0): kx = 1
    if (ky == 0): ky = 1
    if (kz == 0): kz = 1
    if(species in high_k_list):
        _positions_line += 'kpoint_mesh = "23 23 23"\n'
    else:
        _positions_line += 'kpoint_mesh = "%d %d %d"\n'%(kx, ky, kz)

#    _positions_line += 'wavefunction_grid = "%d %d %d"\n'%(nx, ny, nz)
    _positions_line += grids[species] + '"'


    
    jobfile = open('job.bat', 'w') 
    jobfile.write('rm */*.log */*.options\n')
    #jobrun = 'mpirun -n 24 --bind-to core ~/bin/rmg-cpu input\n'
    jobrun = 'srun -AMAT189_crusher --ntasks=256 -u --gpus-per-node=8 --ntasks-per-node=64 --cpu-bind=sockets  ../../rmg-cpu-crusher input\n'

    epps = '\n'
    
    pp_type = '\n'
    if(species in pseudo_uspp):
        pp_type = 'internal_pseudo_type="ultrasoft"\n'
        pp_type += 'potential_grid_refinement="3"\n'

    if(species in pseudo_nc):
        pp_type = 'internal_pseudo_type="nc_accuracy"\n'

    if(species in pseudo_sl):
        pp_type = 'use_bessel_projectors="true"\n'

    mg_line = '\n'
    if(species in run_mg):
        mg_line += 'kohn_sham_solver="multigrid"\n'
        mg_line += 'charge_density_mixing="0.3"\n'
        mg_line += 'charge_mixing_type="Linear"\n'
        mg_line += 'potential_acceleration_constant_step="1.5"\n'
    else:
        mg_line += 'kohn_sham_solver="davidson"\n'
        mg_line += 'charge_density_mixing="0.15"\n'
        mg_line += 'charge_mixing_type="Broyden"\n'
        mg_line += 'potential_acceleration_constant_step="0.0"\n'

    if(species in high_cutoff):
        mg_line += 'kohn_sham_fd_order="10"\n'

    temp_line = 'occupation_electron_temperature_eV = "1.0e-2"\n'
    if(species in high_temp):
        temp_line = 'occupation_electron_temperature_eV = "1.0e-1"\n'

#    if(species in fft_fd):
#        mg_line += 'kohn_sham_ke_fft="true"\n'

    pp_line = '\n'
    if species in external_pps:
        pp_line += 'pseudopotential="' + species + ' ' + external_pps[species] + '"\n'
        pp_line += 'pseudo_dir="../../pseudo"\n'

    extra_line = '\n'
    if species in extra_opts:
        extra_line = extra_opts[species] + '"\n'

    veca = [0.0,0.0,0.0]
    vecb = [0.0,0.0,0.0]
    vecc = [0.0,0.0,0.0]

    for j in range(3):
        veca[j] = crmg.cell.latticevectors[0][j] * crmg.cell.lengthscale
        vecb[j] = crmg.cell.latticevectors[1][j] * crmg.cell.lengthscale
        vecc[j] = crmg.cell.latticevectors[2][j] * crmg.cell.lengthscale
    for vol in volume_lists:
        dir_name = 'volume_' + str(vol)
        jobfile.write('cd ' + dir_name + '\n')
        jobfile.write(jobrun)
        jobfile.write('cd ..\n')
        if not os.path.exists(dir_name): os.mkdir(dir_name)
        f = open(dir_name+'/input', 'w')
        cvol = pow(vol, 1.0/3.0);
        lattice_line = 'lattice_vector = "\n%16.12f %16.12f %16.12f\n%16.12f %16.12f %16.12f\n%16.12f %16.12f %16.12f"\n'%(veca[0]*cvol,veca[1]*cvol,veca[2]*cvol,vecb[0]*cvol,vecb[1]*cvol,vecb[2]*cvol,vecc[0]*cvol,vecc[1]*cvol,vecc[2]*cvol)
        f.write(lattice_line)
#        f.write('a_length="%16.8f"\n'%(a0*vol))
#        f.write('b_length="%16.8f"\n'%(b0*vol))
#        f.write('c_length="%16.8f"\n'%(c0*vol))
        f.write(pp_line)
        f.write(extra_line)
        f.write(temp_line)
        f.write(_positions_line)
        f.write(default_line);
        f.write(pp_type)
        f.write(mg_line)
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

    nx = int ((nx + 1)/4) * 4
    ny = int((ny + 1)/4) * 4
    nz = int((nz + 1)/4) * 4

    
    spacing_list =[a0/nx, b0/ny, c0/nz]
    max_sp = max(spacing_list)
    min_sp = min(spacing_list)
    anisotropy = max_sp/min_sp
    max_iter = 0
    while(anisotropy >1.1 and max_iter < 10):
        max_iter +=1
        nx += 4
        xspacing = a0/nx
        ny = b0/xspacing
        ny = int((ny + 1)/4) * 4
        nz = c0/xspacing
        nz = int((nz + 1)/4) * 4
        spacing_list =[a0/nx, b0/ny, c0/nz]
        max_sp = max(spacing_list)
        min_sp = min(spacing_list)
        anisotropy = max_sp/min_sp
        print (anisotropy)
        

    return(nx, ny, nz)

def read_cif(filename, cif_blocks):

    f = open(filename, 'r')
    all_lines = f.readlines()
    f.close
    
    # process blocks from cif file and stored the values in a dictionay 

    
    num_loops = 0
    for line in all_lines:
        if('loop' in line): num_loops +=1

    loop_lines =[[] for i in range(num_loops)]

    inside_loop = 0;
    for line in all_lines:
        if(inside_loop > 0 and not ('loop' in line) and len(line.strip()) >0 ): 
            loop_lines[inside_loop-1].append(line.strip())
        if('loop' in line): inside_loop +=1

    for line in loop_lines[0]:
        if(line[0] == '_'): dic_item = line
        else: 
           cif_blocks[dic_item] = line
        
    list_keys = []
    for line in loop_lines[1]:
        if(line[0] == '_'): 
            list_keys.append(line)
    values_list = [[] for i in range(len(list_keys))]
    for line in loop_lines[1]:
        if(line[0] != '_'): 
            for i in range(len(list_keys)):
                values_list[i].append(line.split()[i])

    for i in range(len(list_keys)):
        cif_blocks[list_keys[i]] = values_list[i]
        
    for line in all_lines:
        tem = line.strip()
        if(len(tem) == 0): continue
        if(tem[0] == '_'):
           tem1 = tem.split
           cif_blocks[tem.split()[0]] = tem.split()[1]
        if('loop' in tem): break  


if __name__ == '__main__':
    main()
