#!/usr/bin/python3

import sys
import os
import math
import numpy
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

energy_range = 10.0
delta_e = 0.001
gauss_sigma = 50

if len(sys.argv) >1:
    f = open(sys.argv[1], 'r')
else:
    print("use: tddft_post.py dipole.dat spectra.dat")
    print("input file: dipole.dat or *current.dat: dipole(current) moments vs time from RMG")
    print("output file: spectra.dat, the absorption intensity")
    exit()

if len(sys.argv) >2:
    fout = open(sys.argv[2], 'w')
else:
    fout = open('spectra.dat', 'w')
    
#one_direction = 1, only count the dipole of efield direction 
#              = 0: average 3 directions
one_direction = 1


# maximum steps to be postprocesing
n_data_max =20000
#interpolation of energy points

all_lines = f.readlines()
f.close

line = all_lines[1].split()
efield_x = abs(float(line[5]))
efield_y = abs(float(line[6]))
efield_z = abs(float(line[7]))
print("efiled", efield_x, efield_y, efield_z)
if efield_x > 0.0: 
    if efield_y > 0.0: print ('warning: both x and y field >0.0')
    if efield_z > 0.0: print ('warning: both x and z field >0.0')
    col_pick = 1
elif efield_y > 0.0: 
    if efield_z > 0.0: print ('warning: both y and z field >0.0')
    col_pick = 2
elif efield_z > 0.0:
    col_pick = 3
else:
    print ('all x, y, z field == 0')

efield = max(efield_x, efield_y, efield_z)

dt = float(all_lines[4].split()[0])-float(all_lines[3].split()[0])
period = 2.0 * 3.1415926/delta_e * 27.2114

n_data_tot = int(period/dt)
print(n_data_tot)
n_data = min(len(all_lines)-3, n_data_max)
dipole = numpy.zeros(n_data_tot)

for i in range(n_data):
    line = all_lines[i + 3].split()
    if(one_direction):
        dipole[i] = float(line[col_pick]) 
    else:
        dipole[i] = float(line[1]) +float(line[2]) + float(line[3])

      
mean_dipole =  numpy.mean(dipole[0:n_data])

for i in range(n_data):
    dipole[i] = (dipole[i] - mean_dipole) 
#for i in range(n_data):
#    dipole[i] = dipole[i] * math.exp(-float(i)/float(n_data) *5.0)
    

fw = numpy.fft.rfft(dipole, n_data_tot)
  
num_freq = n_data_tot//2 +1

wmin = 0.0
wmax = num_freq * delta_e
x = numpy.linspace(wmin, wmax, num_freq)



spectrum_0 = numpy.zeros(len(fw))
spectrum_2 = numpy.zeros(len(fw))
#  rotate spectrum as in NWCHEM
num_epoint = 0
for i in range(len(fw)):
    re = fw[i].real
    im = fw[i].imag
    r = math.sqrt (re**2 + im**2)
    angle = abs (math.atan2 (im, re))
    if angle > 3.1415927:
      print(angle)
      raise Exception ("atan2 out of range")
    fw[i] = r * abs(math.cos(angle)) + abs(r*math.sin(angle)) * 1j
    if x[i] < energy_range:
        spectrum_0[i] = fw[i].imag/efield
        spectrum_2[i] = fw[i].real/efield
        num_epoint +=1
spectrum_0.resize(num_epoint)
spectrum_2.resize(num_epoint)
spectrum_1 = gaussian_filter1d(spectrum_0, gauss_sigma)
spectrum_3 = gaussian_filter1d(spectrum_2, gauss_sigma)

if "current" in all_lines[2]:
    for i in range(1,num_epoint):
        fout.write("%f  %f\n"%(i*delta_e, spectrum_1[i]/float(i*delta_e)))
    #fout.write("&")
    #for i in range(1,num_epoint):
    #    if i* delta_e >2.0:
    #        fout.write("%f  %f\n"%(i*delta_e, spectrum_3[i]/float(i*delta_e)))
if "BeryPhase" in all_lines[2]:
    for i in range(1,num_epoint):
        fout.write("%f  %f\n"%(i*delta_e, spectrum_1[i]))
    fout.write("&")
    for i in range(1,num_epoint):
        fout.write("%f  %f\n"%(i*delta_e, spectrum_3[i]))
elif "dipole" in all_lines[2]:
    for i in range(1,num_epoint):
        fout.write("%f  %f\n"%(i*delta_e, spectrum_1[i]))
    #fout.write("&")
    #for i in range(1,num_epoint):
    #    fout.write("%f  %f\n"%(i*delta_e, spectrum_3[i]))
