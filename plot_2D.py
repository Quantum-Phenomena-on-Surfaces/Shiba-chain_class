#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 15:01:48 2022

@author: cristina
"""

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
plt.rcParams['figure.figsize'] = [24, 5]
import numpy as np

def plot_2D(z, z_z, z_x, energy, e, e_z, e_x, vv, N_x):
    
    eV = 27.2114
    
    fig = plt.figure(2)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    im1 = ax1.imshow(z/eV, aspect= 'equal', cmap = 'viridis')
    ax1.title.set_text('PDOS at E = %1.4f meV' %energy)
    fig.colorbar(im1, ax=ax1, orientation='horizontal')
    
    im2 = ax2.imshow(z_z/eV, aspect= 'equal', cmap = 'RdBu_r')
    ax2.title.set_text(r'$\rho_{\uparrow} - \rho_{\downarrow}$')
    fig.colorbar(im2, ax=ax2, orientation='horizontal')
    
    im3 = ax3.imshow(z_x/eV, aspect= 'equal', cmap = 'RdBu_r')
    ax3.title.set_text(r'$\rho_{x}$')
    fig.colorbar(im3, ax=ax3, orientation='horizontal')
    
    #save Fig   
    plt.savefig('2D_map_z.png')
    
    #savedata
    np.savetxt('z_total.txt', z/eV)
    np.savetxt('z_z.txt', z_z/eV)
    np.savetxt('z_x.txt', z_x/eV)
    
    fig2 = plt.figure(3)
    ax1 = fig2.add_subplot(131)
    ax2 = fig2.add_subplot(132)
    ax3 = fig2.add_subplot(133)
    
    im1 = ax1.imshow(e/eV, aspect='auto', cmap = 'viridis', extent=[vv[0], vv[-1], 1, N_x]) 
    cb1 = fig.colorbar(im1, ax=ax1, orientation='horizontal')
    ax1.set_xlabel('Energy (meV)')
    cb1.set_label('PDOS (1/eV)')
    
    
    im2 = ax2.imshow(e_z/eV, aspect='auto', cmap = 'viridis', extent=[vv[0], vv[-1], 1, N_x]) 
    cb2 = fig.colorbar(im2, ax=ax2, orientation='horizontal')
    ax2.set_xlabel('Energy (meV)')
    cb2.set_label(r'$\rho_{\uparrow} - \rho_{\downarrow}$')
    
    im3 = ax3.imshow(e_x/eV, aspect='auto', cmap = 'viridis', extent=[vv[0], vv[-1], 1, N_x]) 
    cb3 = fig.colorbar(im3, ax=ax3, orientation='horizontal')
    ax3.set_xlabel('Energy (meV)')
    cb3.set_label(r'$\rho_{x}$')
    
    #save Fig   
    plt.savefig('2D_map_zepng')
    
    #savedata
    np.savetxt('e_total.txt', e/eV)
    np.savetxt('e_z.txt', e_z/eV)
    np.savetxt('e_x.txt', e_x/eV)    