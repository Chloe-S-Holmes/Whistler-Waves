#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:45:11 2024

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
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FuncFormatter

"""
***
***
DONT FORGET TO UPDATE PRESSURE AND FREQUENCY WHEN YOU CHANGE FILE
***
***
"""
pressure='1p5'
frequency='100'
fname = '/home/chloe/Downloads/JuniorYear/PlasmaLab/Experiments/Lab3/data/Lab3/Bdot_n10_n60_100mhz.hdf5'
def csv_output_file(string, k):
    if k == True:
        return f'/home/chloe/Downloads/JuniorYear/PlasmaLab/Experiments/Lab3/data/extracted_parameters/{pressure}mTorr/k_values/{frequency}Hz.csv'
    else:
        return f'/home/chloe/Downloads/JuniorYear/PlasmaLab/Experiments/Lab3/data/extracted_parameters/{pressure}mTorr/fitting_params/{string}{frequency}Hz.csv'

data_set = DataSet(fname)
data_set.print_tree()

Bx = data_set.read_data(1)
By = data_set.read_data(2)
time = data_set.read_time()
pos = data_set.read_positions()
z_pos = -np.array([x[1] for x in pos])
#time = time[100:700]
#Bx=Bx[:,100:700]

def plot_overlapping_B(index, axis='x'):
    plt.title(f'Magnetic field vs time for {-pos[-(index+1)][1]}cm and {-pos[index][1]}cm')
    plt.xlabel('index of time array')
    
    time_index = [i for i in range(len(time))]
    if axis =='x':
        plt.plot(time_index, Bx[index])
        plt.plot(time_index, Bx[-(index+1)])
        plt.ylabel('Bx (Arb)')
    elif axis == 'y':
        plt.plot(time_index, By[index])
        plt.plot(time_index, By[-(index+1)])
        plt.ylabel('By (Arb)')
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(n=5))
    plt.grid(axis='x', which='both')
    plt.show()

def bandpass_filter(signal, f_cutoff, sample_freq=2.5e9, order=2):
    sos = butter(order, f_cutoff, 'bp', fs=sample_freq, output='sos')
    filtered_signal = sosfiltfilt(sos, signal)
    return filtered_signal
n_pos=3
plot_overlapping_B(n_pos)
n_time=[215,220,225,230,235,240]

Bx_filt = bandpass_filter(Bx, [1e6, 140e6])


Bx_truncated=Bx_filt[n_pos:-(n_pos+1)]
z_pos_truncated=z_pos[n_pos:-(n_pos+1)]

def damped_waveform(z, A, ki, kr, phi):
    return A*np.exp(-ki*z)*np.cos(kr*z+phi)

def plot_waveform_vs_pos(params, z, B, n_time):
    plt.xlabel('Distance from whistler wave exciter (cm)')
    plt.ylabel('Bx (Arb)')
    t = '{:.2g}'.format(time[n_time])
    plt.title(f'Position Variance of the Whistler Wave at Time {t}s')
    z_fit=np.linspace(z[0],z[-1],100)
    B_fit = [damped_waveform(i,*params) for i in z_fit]
    plt.plot(z_fit,B_fit, label='curve fit', color='orange')
    plt.plot(z, B, label='B data', color='blue')
    plt.legend()
    plt.show()



bounds=([1e-5, 1e-5, 1e-6, -np.pi/2],[10, 10, 10, np.pi/2])
p0=[(3, 4e-1, 6e-2, 0),
    (3, 1e-1, 1e-5, 0),
    (3, 3e-1, 3e-1, 0),
    (3, 1e-1, 8e-1, 0),
    (3, 1e-1, 5e-1, 0),
    (3, 2e-3, 5e-1, 0)]

i=5
B_of_z=Bx_truncated[:,n_time[i]]
z=z_pos_truncated
popt, cov = optimize.curve_fit(damped_waveform, z, B_of_z, bounds=bounds, p0=p0[i])
plot_waveform_vs_pos(popt, z, B_of_z, n_time[i])
save_params=False


if len(p0) != len(n_time):
    print('p0 is wrong length')

def find_all_k():
    kr=[]*len(n_time)
    ki=[]*len(n_time)
    for i in range(len(n_time)):
        B_of_z=Bx_truncated[:,n_time[i]]
        z=z_pos_truncated
        popt, cov = optimize.curve_fit(damped_waveform, z, B_of_z, bounds=bounds, p0=p0[i])
        kr.append(popt[2])
        ki.append(popt[1])
    return kr, ki

if save_params==True:
    kr,ki=find_all_k()
    k_values = [kr,ki]
    np.savetxt(csv_output_file('kvalues',True), k_values, delimiter=',', fmt='%s')
    np.savetxt(csv_output_file('n_time',False), n_time, delimiter=',', fmt='%s')
    np.savetxt(csv_output_file('bounds',False), bounds, delimiter=',', fmt='%s')
    np.savetxt(csv_output_file('p0',False), p0, delimiter=',', fmt='%s')
    print(k_values)
    





















