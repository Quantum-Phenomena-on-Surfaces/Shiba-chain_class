#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 14:29:41 2022

@author: cristina
"""

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
plt.rcParams['figure.figsize'] = [8, 12]
import numpy as np

def plot_spectra(vv, spect_atom, spect_atom_up, spect_atom_dn, profile, energy):       
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.plot(vv, spect_atom, label = 'total')
    ax1.plot(vv, spect_atom_up, label = 'up')
    ax1.plot(vv, spect_atom_dn, label = 'down')
    
    ax1.set_ylabel('PDOS')
    ax1.set_xlabel('Energy (meV)')
    ax1.title.set_text('First Magnetic atom')
    ax1.set_xlim(min(vv), max(vv))
    ax1.legend()
    
    
    ax2.plot(profile, '--o')
    ax2.set_ylabel('Profile at E = %f meV along x' %energy)
    ax2.set_xlabel('Atomic site')
    
    #save data
    
    #plt.savefig('spectrum.png')
    