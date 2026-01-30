#!/usr/bin/env python3
import numpy as np
import sys
import os.path


if __name__ == '__main__':
    if os.path.exists('input.ref_spin0_current.dat'):
        ref_name = 'input.ref_spin0_current.dat'
        new_name = 'input.00_spin0_current.dat'
    elif os.path.exists('input.ref_spin0_dipole.dat'):
        ref_name = 'input.ref_spin0_dipole.dat'
        new_name = 'input.00_spin0_dipole.dat'
  
    ref_data = np.loadtxt(ref_name, skiprows = 3)
    new_data = np.loadtxt(new_name, skiprows = 3)

    if(ref_data.shape != new_data.shape):
        print(ref_data.shape, new_data.shape)
    diff = ref_data[:,1:4] - new_data[:,1:4]
    abs_sum_diff = np.sum(np.abs(diff))
    abs_sum_ref = np.sum(np.abs(ref_data[:, 1:4]))
    current_pass = abs_sum_diff/abs_sum_ref < 0.02

    ref_name = 'Epsilon/absorb_spin0.dat_ref'
    new_name = 'Epsilon/absorb_spin0.dat'
    
  
    ref_data = np.loadtxt(ref_name, skiprows = 1)
    new_data = np.loadtxt(new_name, skiprows = 1)

    if(ref_data.shape != new_data.shape):
        print(ref_data.shape, new_data.shape)
    diff = ref_data[:,1:4] - new_data[:,1:4]
    abs_sum_diff = np.sum(np.abs(diff))
    abs_sum_ref = np.sum(np.abs(ref_data[:, 1:4]))
    absorb_pass = abs_sum_diff/abs_sum_ref < 0.01

    ref_name = 'Epsilon/epsilon_spin0.dat_ref'
    new_name = 'Epsilon/epsilon_spin0.dat_ref'
    
  
    ref_data = np.loadtxt(ref_name, skiprows = 1, comments = "&")
    new_data = np.loadtxt(new_name, skiprows = 1, comments = "&")

    if(ref_data.shape != new_data.shape):
        print(ref_data.shape, new_data.shape)
    diff = ref_data[:,1:4] - new_data[:,1:4]
    abs_sum_diff = np.sum(np.abs(diff))
    abs_sum_ref = np.sum(np.abs(ref_data[:, 1:4]))
    epsilon_pass = abs_sum_diff/abs_sum_ref < 0.01

    if(current_pass and absorb_pass and epsilon_pass):
        print("pass")
    else:
        if not current_pass:
            print("current failed")
        if not absorb_pass:
            print("absorb failed")
        if not epsilon_pass:
            print("epsilon failed")
