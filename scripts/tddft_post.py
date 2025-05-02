#!/usr/bin/python3

import sys
import os
import math
import numpy
import matplotlib.pyplot as plt

if len(sys.argv) >1:
    f = open(sys.argv[1], 'r')
else:
    print("use: tddft_post.py dipole.dat spectra.dat")
    print("input file: dipole.dat: dipole moments vs time from RMG")
    print("output file: spectra.dat, the absorption intensity")
    exit()

if len(sys.argv) >2:
    fout = open(sys.argv[2], 'w')
else:
    fout = open('spectra.dat', 'w')
    
energy_range = 10.0
#one_direction = 1, only count the dipole of efield direction 
#              = 0: average 3 directions
one_direction = 1


# maximum steps to be postprocesing
n_data_max = 5000
#interpolation of energy points
ndup = 4

# nsmooth = 2: no smooth
nsmooth = 2

all_lines = f.readlines()
f.close

line = all_lines[1].split()
efield_x = abs(float(line[2]))
efield_y = abs(float(line[3]))
efield_z = abs(float(line[4]))
  
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
  
      
n_data = min(len(all_lines) -3, n_data_max)

dipole = numpy.zeros(n_data)
t = numpy.zeros(n_data)

for i in range(n_data):
    line = all_lines[i + 3].split()
    t[i] = float(line[0])
    if(one_direction):
        dipole[i] = float(line[col_pick]) 
    else:
        dipole[i] = float(line[1]) +float(line[2]) + float(line[3])

      
mean_dipole =  numpy.mean(dipole)

for i in range(n_data):
    dipole[i] = (dipole[i] - mean_dipole) 
    

fw = numpy.fft.rfft(dipole, n_data)
  
num_freq = n_data//2 +1
dt = t[1] - t[0]
period = n_data * dt
dw = 2.0 * 3.1415926/period * 27.2114

wmin = 0.0
wmax = num_freq * dw
x = numpy.linspace(wmin, wmax, num_freq)



spectrum_0 = numpy.zeros(len(fw))
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
    fw[i] = r * math.cos(angle) + r*math.sin(angle) * 1j
    if x[i] < energy_range:
        spectrum_0[i] = fw[i].imag/efield
        num_epoint +=1
spectrum_0.resize(num_epoint)
fw1 = numpy.fft.fft(spectrum_0, num_epoint)
fw2 = numpy.zeros(num_epoint *ndup) + 0.0* 1j
fw2[0] = fw1[0]
for i in range(1, len(fw1)//nsmooth+1):
    fw2[i] = fw1[i]
    fw2[num_epoint*ndup-i] = fw1[num_epoint -i]
fw = numpy.fft.ifft(fw2, num_epoint * ndup)
for i in range(len(fw)):
    fout.write("%f  %f\n"%(i*dw/ndup, (ndup * fw[i].real/n_data)))
fout.write("&\n")
