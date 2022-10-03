#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 14:29:41 2022

@author: cristina
"""

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
plt.rcParams['figure.figsize'] = [8, 10]
import numpy as np

def plot_spectra(vv, spect_atom, spect_atom_up, spect_atom_dn, profile, energy):       
    
    eV = 27.2114
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.plot(vv, spect_atom/eV, label = 'total')
    ax1.plot(vv, spect_atom_up/eV, label = 'up')
    ax1.plot(vv, spect_atom_dn/eV, label = 'down')
    
    ax1.set_ylabel('PDOS (1/eV)')
    ax1.set_xlabel('Energy (meV)')
    ax1.title.set_text('First Magnetic atom')
    ax1.set_xlim(min(vv), max(vv))
    ax1.legend()
    
    
    ax2.plot(profile/eV, '--o')
    ax2.set_ylabel('Profile at E = %f meV along x' %energy)
    ax2.set_xlabel('Atomic site')
    
    #save Fig   
    plt.savefig('spectra_profile.png')
    
    #save data
    np.savetxt('spect_total.txt', spect_atom/eV)
    np.savetxt('spect_up.txt', spect_atom_up/eV)
    np.savetxt('spect_down.txt', spect_atom_dn/eV)
    np.savetxt('vv.txt', vv)
    
    