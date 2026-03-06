#!/usr/bin/python3

import sys
import os
import math
import numpy
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

col_pick = [1,1,0] 
dumping = 0.0
energy_range = [0.0,20.0]
delta_e = 0.001
gauss_sigma = 40
# maximum steps to be postprocesing
n_data_max = 80000

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


#interpolation of energy points

all_lines = f.readlines()
f.close

line = all_lines[1].split()
efield = [0.0, 0.0, 0.0]
efield[0] =  abs(float(line[5]))
efield[1] =  abs(float(line[6]))
efield[2] =  abs(float(line[7]))
efield_max = max(efield)
print("efiled", efield,efield_max)

dt = float(all_lines[4].split()[0])-float(all_lines[3].split()[0])
period = 2.0 * 3.1415926/delta_e * 27.2114

n_data_tot = int(period/dt)
print(n_data_tot)

n_data = min(len(all_lines)-3, n_data_max)

num_freq = n_data_tot//2 +1
spectrum_i = numpy.zeros(num_freq)
spectrum_r = numpy.zeros(num_freq)
for ixyz in range(3):
    if col_pick[ixyz] == 0: continue 
    dipole = numpy.zeros(n_data_tot)

    for i in range(n_data):
        line = all_lines[i + 3].split()
        dipole[i] = float(line[ixyz+1]) 
          
    mean_dipole =  numpy.mean(dipole[0:n_data])

    for i in range(n_data):
        dipole[i] = (dipole[i] - mean_dipole) 
    for i in range(n_data):
        dipole[i] = dipole[i] * math.exp(-float(i)/float(n_data) *dumping)
        
    if "current" in all_lines[2]:
        dipole1 = dipole[0:n_data].cumsum()
        x = numpy.linspace(0, 1, n_data)
        linfit = numpy.poly1d(numpy.polyfit(x, dipole1, 1))
        dipole[0:n_data] = dipole1 - linfit(x)

    fw = numpy.fft.rfft(dipole, n_data_tot)
      
    wmin = 0.0
    wmax = num_freq * delta_e
    x = numpy.linspace(wmin, wmax, num_freq)



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
        if x[i] <= energy_range[1]:
            spectrum_i[i] += fw[i].imag/efield_max
            spectrum_r[i] += fw[i].real/efield_max
            num_epoint +=1

spectrum_r.resize(num_epoint)
    
spectrum_rs = gaussian_filter1d(spectrum_r, gauss_sigma)
for i in range(1,num_epoint):
    if i*delta_e >= energy_range[0]:
        fout.write("%f  %f\n"%(i*delta_e, spectrum_rs[i]))
