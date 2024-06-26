#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 23:08:56 2024

@author: chloe
"""
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from scipy import optimize

frequencies = [40,45,50,55,60,65,70,75,80,85,90,96,98,100,102,104,106,108,110,112,114,116,118,120]

kr=[]
ki=[]
err_kr=[]
err_ki=[]

for i in range(len(frequencies)):
    filename = f'/home/chloe/Downloads/JuniorYear/PlasmaLab/Experiments/Lab3/data/extracted_parameters/1p5mTorr/k_values/{frequencies[i]}Hz.csv'
    data = pd.read_csv(filename, header = None)
    
    kr_values = np.array(data.iloc[0])
    ki_values = np.array(data.iloc[1])
    
    kr_mean = np.mean(kr_values, axis=0)
    kr_err = np.std(kr_values, axis=0)
    ki_mean = np.mean(ki_values, axis=0)
    ki_err = np.std(ki_values, axis=0)
    
    kr.append(kr_mean)
    err_kr.append(kr_err)
    ki.append(ki_mean)
    err_ki.append(ki_err)

def dispersion_relation(w,de,wce):
    kr = (1/de)*(w/(np.abs(wce)-w))**(1/2)
    return kr



w = 10**6*2*np.pi*np.array(frequencies)

popt, cov = optimize.curve_fit(dispersion_relation, w,kr, bounds=([1e-1, 1e6],[100,1e10]), p0=[2.42, 1e9])
std = np.sqrt(np.diag(cov))
print('extracted de,wce')
print(popt,std)

plt.errorbar(kr, w, xerr=err_kr, yerr=2*np.pi*1, fmt='none',capsize=5)
plt.scatter(kr, w, s= 10)
plt.plot(dispersion_relation(w, *popt), w, label = 'fitted dispersion relation')

ne=8.7e10
Te=3.5

B=75.6
de=5.31e5*(ne**(-1/2))
wce=1.76e7*B

print('predicted de, wce')
print(de,wce)
plt.title('Dispersion Relation of Whistler Waves at 353V 1.5mTorr')
plt.ylabel(f'$\omega$ (rad/s)')
plt.xlabel(f'$k_r$ $(cm^{-1})$')
plt.plot(dispersion_relation(w,de,wce), w, label='predicted dispersion relation')
plt.legend(loc='upper left')
plt.show()


de=popt[0]
wce=popt[1]

def imaginary_relation(w,v):
    ki = (v/(2*de*(wce-2)))*((w*(wce-w))/((wce-w)**2-v**2))**(1/2)
    return ki

p=0.2
kb=1.38e-23
T=295
kbT=kb*T
sigma=1e-20

de=5.31e5*(ne**(-1/2))
wce=1.76e7*B

vth=4.19e7*(Te**(1/2))
v=vth/mfp
print('predicted v (/e6)')
print(v/10e6)

nn=p/(kbT)
mfp=1/(nn*sigma)

popt, cov = optimize.curve_fit(imaginary_relation, w,ki, bounds=(1e5,1e9), p0=3.56e6)
std = np.sqrt(np.diag(cov))
print('extracted v (/e6)')
print(popt/10e6,std/10e6)

plt.errorbar(ki, w, xerr=err_kr, yerr=2*np.pi*1, fmt='none',capsize=5)
plt.scatter(ki, w, s= 10)
plt.xlabel(f'$k_i$ $cm^{-1}$')
plt.ylabel(f'$\omega$ (rad/s)')
plt.title('Imaginary Component of the Wavenumber vs Frequency')
plt.plot(imaginary_relation(w, *popt), w, label = 'fitted k{Im}')
plt.plot(imaginary_relation(w,v),w, label='predicted k{Im}')
plt.legend()
plt.show()







