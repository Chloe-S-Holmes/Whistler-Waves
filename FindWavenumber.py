#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:04:36 2024

@author: chloe
"""

import sys
sys.path.append('/home/chloe/Downloads/JuniorYear/PlasmaLab/Experiments/Lab3/code')
from lecroydaq import DataSet
import matplotlib.pylab as plt
ax = plt.gca()
from scipy.signal import savgol_filter, find_peaks, butter, sosfiltfilt
import numpy as np
from scipy import optimize

fname = '/home/chloe/Downloads/JuniorYear/PlasmaLab/Experiments/Lab3/data/Lab3/Bdot_n10_n60_120mhz.hdf5'
data_set = DataSet(fname)
data_set.print_tree()

Bx = data_set.read_data(1)
time = data_set.read_time()
pos = data_set.read_positions()

def bandpass_filter(signal, f_cutoff, sample_freq=2.5e9, order=2):
    sos = butter(order, f_cutoff, 'bp', fs=sample_freq, output='sos')
    filtered_signal = sosfiltfilt(sos, signal)
    return filtered_signal

def B_of_z(z, A, ki, kr, phi):
    return A*np.exp(-ki*z)*np.cos(kr*z+phi)

pos_max=40
time_sample=200

plt.plot(Bx[pos_max])
plt.plot(Bx[0])
plt.show()


Bx_filt = bandpass_filter(Bx, [1e6, 140e6])
posz = [pos[i][1] for i in range(51)]
Bx_time = Bx_filt[:, 200]
z = np.abs(posz[:])
plt.plot(z, Bx_time)

popt, cov = optimize.curve_fit(B_of_z, z, Bx_time, bounds=([1e-5, 1e-4, 1e-1, -np.pi/2],[10, 1, 100, np.pi/2]))
print(popt)
#plt.plot(z, 2*np.exp(-1e-3*z))
plt.plot(z, B_of_z(z, *popt))

kr = popt[2]