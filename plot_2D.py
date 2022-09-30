#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 15:01:48 2022

@author: cristina
"""

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
plt.rcParams['figure.figsize'] = [24, 4]
import numpy as np

def plot_2D(z, z_z, z_x, energy, e, e_z, e_x, vv, N_x):
    
    fig = plt.figure(2)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    im1 = ax1.imshow(z, aspect= 'equal', cmap = 'viridis')
    ax1.title.set_text('PDOS at E = %1.4f meV' %energy)
    fig.colorbar(im1, ax=ax1)
    
    im2 = ax2.imshow(z_z, aspect= 'equal', cmap = 'RdBu_r')
    ax2.title.set_text(r'$\rho_{\uparrow} - \rho_{\downarrow}$')
    fig.colorbar(im2, ax=ax2)
    
    im3 = ax3.imshow(z_x, aspect= 'equal', cmap = 'RdBu_r')
    ax3.title.set_text(r'$\rho_{x}$')
    fig.colorbar(im3, ax=ax3)
    
    
    fig2 = plt.figure(3)
    ax1 = fig2.add_subplot(131)
    ax2 = fig2.add_subplot(132)
    ax3 = fig2.add_subplot(133)
    
    im1 = ax1.imshow(e, aspect='auto', cmap = 'viridis', extent=[vv[0], vv[-1], 1, N_x]) 
    fig.colorbar(im1, ax=ax1, orientation='horizontal')
    
    im2 = ax2.imshow(e_z, aspect='auto', cmap = 'viridis', extent=[vv[0], vv[-1], 1, N_x]) 
    fig.colorbar(im2, ax=ax2, orientation='horizontal')
    
    im3 = ax3.imshow(e_x, aspect='auto', cmap = 'viridis', extent=[vv[0], vv[-1], 1, N_x]) 
    fig.colorbar(im3, ax=ax3, orientation='horizontal')
    
    
    
    