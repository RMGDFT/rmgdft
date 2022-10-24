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

pi = 3.141592653589793238

# The base value is used for ncpp
grid_spacing_base = 0.12 # in unit of Angstrom
volume_lists = [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06]
pseudo_extra = ["Ar", "Bi", "He", "Kr", "Lu", "Ne", "Po", "Rn", "Xe"]

# Rmg default pseudo is NCPP. If a species is in this list then USPP will be used.
pseudo_uspp = ["Ag","O","Mn","Ga","Mg","Cr","Mo","Tc","Zn","F","Cl","Ni","Sc","Y","B","Co","Po","Ge","Au","As","Si","V","Os","I","Li","Cs","In","K","Ru","Ti","Hg","Tl","H","Sb","Se","S","Rh","W","Re","Ir","Zr","Al","Cd","Br","Be","Ba","P","Nb","Na","Pt","Rb","Pb","Pd","Sr","Ta","Sn","Ca","C","Cu","Hf","Te"];

# Default solver is davidson but multigrid will be used if a species is in this list
run_mg = ["O"];

FM_list = ["Fe", "Co", "Ni"]
AFM_list = ["Cr", "Mn", "O"]


k_delta = 0.08  # in unit of (Anstrom^-1)
k_parall = 4;
pp = 'MIXED'
#pp = 'ONCV'

def main():
    jobline = """
#!/bin/bash
cp ~/bin/rmg-cpu .
"""

    jobline_desktop = ''
#    line_run = 'aprun -n 128 -N 8 -d 4 -cc numa_node ../../rmg-cpu input&\n'
    line_run = 'srun -AMAT189_crusher --ntasks=32 -u -c8 --gpus-per-node=8 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=single:1 --cpu-bind=sockets  ../../rmg-gpu-crusher input\n'

    
    shell_script_line=''
    dirname = 'Delta_benchmark_dk_%4.2f_hx_%4.2f_%s'%(k_delta, grid_spacing_base, pp)
    if not os.path.exists(dirname): os.mkdir(dirname)
    os.chdir(dirname)
    for filename in os.listdir('../../CIFs/'):
        species = filename.split('.')[0]
    
        cif_blocks = {}
        read_cif('../../CIFs/'+filename, cif_blocks)
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

        input_for_rmg(species, cif_blocks, grid_spacing)
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

def input_for_rmg(species, cif_blocks, grid_spacing):



    default_line = """
crds_units = "Angstrom"
lattice_units = "Angstrom"
atomic_coordinate_type = "Cell Relative"
potential_grid_refinement="2"
occupations_type = "MethfesselPaxton"
occupation_electron_temperature_eV = "1.0e-2"
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

    atom_list = []
    for i in range(len(cif_blocks['_atom_site_type_symbol'])):
        atom_name = cif_blocks['_atom_site_type_symbol'][i]
        x = float(cif_blocks['_atom_site_fract_x'][i])
        y = float(cif_blocks['_atom_site_fract_y'][i])
        z = float(cif_blocks['_atom_site_fract_z'][i])
        atom_list.append([atom_name, x, y, z])

    a0 = float(cif_blocks['_cell_length_a'])
    b0 = float(cif_blocks['_cell_length_b'])
    c0 = float(cif_blocks['_cell_length_c'])
    alpha = float(cif_blocks['_cell_angle_alpha'])
    beta = float(cif_blocks['_cell_angle_beta'])
    gamma =float(cif_blocks['_cell_angle_gamma'])
#    if(abs(alpha -60) < 0.1 and abs(beta-60) < 0.1 and abs(gamma-60) <0.1 and abs(a0-b0) < 0.01 and abs(a0-c0) < 0.01):
#        a0 = a0*math.sqrt(2.0)
#        b0 = a0
#        c0 = a0
#    elif(abs(alpha -90) < 0.1 and abs(beta-90) < 0.1 and abs(gamma-120) <0.1 and abs(a0-b0) < 0.01):
#        b0 = b0 * math.sqrt(3.0)
#    elif(abs(alpha -90) < 0.1 and abs(beta-90) < 0.1 and abs(gamma-60) <0.1 and abs(a0-b0) < 0.01):
#        b0 = b0 * math.sqrt(3.0)
    
    veca = [0.0,0.0,0.0]
    vecb = [0.0,0.0,0.0]
    vecc = [0.0,0.0,0.0]
    alpha_r = alpha * pi / 180.0
    beta_r = beta * pi / 180.0
    gamma_r = gamma * pi / 180.0
    veca[0] = a0;
    vecb[0] = b0 * cos(gamma_r)
    vecb[1] = b0 * sin(gamma_r)
    vecc[0] = c0 * cos(beta_r)
    cy = ( cos(alpha_r) - cos(gamma_r) * cos(beta_r) ) / sin(gamma_r)
    vecc[1] = c0 * cy
    vecc[2] = c0*sqrt(1.0-cos(alpha_r)*cos(alpha_r)-cos(beta_r)*cos(beta_r)-cos(gamma_r)*cos(gamma_r)+2*cos(alpha_r)*cos(beta_r)*cos(gamma_r))/sin(gamma_r)

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
    #jobrun = 'mpirun -n 24 --bind-to core ~/bin/rmg-cpu input\n'
    jobrun = 'srun -AMAT189_crusher --ntasks=32 -u -c8 --gpus-per-node=8 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=single:1 --cpu-bind=sockets  ../../rmg-gpu-crusher input\n'

    pp_type = '\n'
    if(species in pseudo_uspp):
        pp_type = 'internal_pseudo_type="ultrasoft"'

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

    pp_line='\n'
#    if(pp == 'ONCV'):
#        pp_line = ' pseudopotential = "%s  ../../../ONCV/%s_ONCV_PBE-1.0.upf"\n'%(atom_list[0][0], atom_list[0][0])
#
#    if(species in pseudo_extra):
#        pp_line ='pseudopotential = "%s ../../../PBE_MT_FHI/%s.pbe-mt_fhi.txt"\n '%(species, species)
#

    for vol in volume_lists:
        dir_name = 'volume_' + str(vol)
        jobfile.write('cd ' + dir_name + '\n')
        jobfile.write(jobrun)
        jobfile.write('cd ..\n')
        if not os.path.exists(dir_name): os.mkdir(dir_name)
        f = open(dir_name+'/input', 'w')
        cvol = pow(vol, 1.0/3.0);
        lattice_line = 'lattice_vector = "\n%f %f %f\n%f %f %f\n%f %f %f"\n'%(veca[0]*cvol,veca[1]*cvol,veca[2]*cvol,vecb[0]*cvol,vecb[1]*cvol,vecb[2]*cvol,vecc[0]*cvol,vecc[1]*cvol,vecc[2]*cvol)
        f.write(lattice_line)
#        f.write('a_length="%16.8f"\n'%(a0*vol))
#        f.write('b_length="%16.8f"\n'%(b0*vol))
#        f.write('c_length="%16.8f"\n'%(c0*vol))
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
        ny = int((ny + 1)/2) * 2
        nz = c0/xspacing
        nz = int((nz + 1)/2) * 2
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
