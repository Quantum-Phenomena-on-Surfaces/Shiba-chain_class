#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 16:11:14 2022

@author: cristina
"""

import numpy as np
import parameters as p     
import Self_Energy_3D as SL
import Free_Green_new as FG
import Dyson as Dy
import detect_peaks as dp
import time

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

class Shiba_chain:    
    def __init__(self):#initialiatise class by calculating total G_chain
        
        t1 = time.time()
           
        pi = np.pi
        d = p.d
        N_atoms = p.N_atoms
        borde = p.borde
        ancho = p.wide
        layers = 1
        lattice_param = p.lattice_param
        a_interatomic = d*lattice_param/0.529

        #load parameters
        alpha = p.alpha
        lamda = (alpha/(2*a_interatomic*0.529))/27.2116
        state = p.state#spin state
        N_period = p.N_period#spin helix 2=AFM
        k_F = p.k_F
        u = p.U#%potential scatt#
        j = p.j#magnetic coupling
        DOS_o = p.DOS
        s = p.s
        Delta = p.delta
        N_omega = p.N_omega
        omega0 = p.omega0
        omega1 = p.omega1
        Dynes = p.Dynes     
        Damping = Dynes/27211.6
        mass_eff=1
        
        '''Asign all parameters'''
        self.a_interatomic = a_interatomic
        self.k_F = k_F
        #self.U = U
        #self.j = J
        self.delta = Delta
        self.dynes = Dynes
        self.rashba = alpha
        self.N_omega = N_omega  
        self.borde = borde
        self.ancho = ancho
        
        N_x = N_atoms + 2*borde#+1 to make the same lattice without impurities
        N_y = ancho    
        
        row = int(N_y/2)
        medio=int(N_x/2)
        
        self.N_x = N_x
        self.N_y = N_y
        self.row = row
        self.medio = medio
        
        '''J vector'''
        J = np.zeros(N_atoms, dtype = float)
        J[::1] = j
        
        '''U vector'''
        U = np.zeros(N_atoms, dtype = float)
        U[::1] = u
        
        self.J = J
        self.U = U
        
        "Spin arragement"
        S = s
        N_periodicty = N_period   
        
        if (state == 'FM'):
            thetaS = np.zeros(N_atoms, dtype = 'float')
            
        elif (state == 'AF'):
            thetaS = np.zeros(N_atoms, dtype = 'float')
            for i in range(int(len(thetaS)/2)):
                thetaS[2*i+1] = pi
            if N_atoms % 2 == 0:
                thetaS[-1] = pi         
        
        elif (state == 'helix'):        
            theta = np.linspace(0.0, 2*pi, N_periodicty+1, dtype = 'float')
            theta = theta[0:N_periodicty]####remove 2pi 
            thetaS = np.zeros(N_atoms, dtype = 'float')
            t = []
        
            if(N_periodicty < N_atoms):
    
                while (len(t)<N_atoms):
                    t = np.concatenate((t,theta))        
                thetaS[0:N_atoms] = t[0:N_atoms]
        
            else:                
                thetaS[0:N_atoms] = theta[0:N_atoms]  
                
        phi = np.zeros(N_atoms)#phi is set to zero
        
        self.thetaS = thetaS
        self.phi = phi

        #########omega vector
        Romega = np.linspace(omega0/27211.6, omega1/27211.6, N_omega)
        vv=Romega*27211.6        
        
        "We calculate the Green's functions and solve Dyson eq"        
        Self = SL.Self_Energy(J, S, thetaS, phi, U, N_atoms, N_x, N_y, layers, borde, lamda)
             
        GG = np.zeros([4 * N_y * N_x*layers , 4 *N_y * N_x*layers, N_omega], dtype=complex)
    
    
        for i_omega in range(N_omega):
        
            omega = Romega[i_omega]    
            #BCS Green's function
            Go = FG.Free_Green(N_x, N_y, layers, omega, Damping, k_F, mass_eff, DOS_o, Delta, a_interatomic)
            
            #Solve Dyson's equation
            gg = Dy.Dyson_eq(Go , Self , N_x, N_y, layers)
        
            GG[:,:, i_omega] = gg
            
        '''Asign resulting matrix'''    
        self.GG = GG
        self.vv = vv
        
        t2 = time.time()
        
        print('The program is finished after', t2 - t1)
        print('Number of atoms', N_y * N_x * layers)
        
        
    def spectra_calc(self):
        pi = np.pi
        N_x = self.N_x
        N_y = self.N_y
        N_omega = self.N_omega
        GG = self.GG
        
        ####total spectra
        spectra = np.zeros([N_y, N_x, N_omega], dtype= 'float')
        spectra_up = np.zeros([N_y, N_x, N_omega], dtype= 'float')
        spectra_down = np.zeros([N_y, N_x, N_omega], dtype= 'float')
        spectra_x =  np.zeros([N_y, N_x, N_omega], dtype= 'float')
        
        for i_atom in range(N_y):
            for j_atom in range(N_x):            
                I = i_atom*N_x + j_atom
        
                for i_omega in range(N_omega):                
                
                    tr = GG[I*4 + 0, I*4 + 0, i_omega] + GG[I*4 + 3, I*4 + 3, - (i_omega+1)]
                    tr_up = GG[I*4 + 0, I*4 + 0, i_omega]
                    tr_dn = GG[I*4 + 3, I*4 + 3, - (i_omega+1)]                
                    trx = GG[I*4 + 2, I*4 + 3, i_omega] + GG[I*4 + 3, I*4 + 2, i_omega]
                
                    spectra[i_atom , j_atom, i_omega]= - (tr.imag)/(pi)
                    spectra_up[i_atom , j_atom, i_omega] = - (tr_up.imag)/(pi)
                    spectra_down[i_atom , j_atom, i_omega] = - (tr_dn.imag)/(pi)
                    spectra_x[i_atom , j_atom, i_omega] = - (trx.imag)/(pi)               
                
                
        '''Asign resulting spectra'''      
        self.spectra = spectra        
        self.spectra_up = spectra_up    
        self.spectra_down = spectra_down    
        self.spectra_x = spectra_x
        
        
    def select_spectra(self, i, j):          
            
        spectrum_atom = self.spectra[i, j, :]#spectrum at the i,j location
        spectrum_atom_up = self.spectra_up[i, j, :]
        spectrum_atom_down = self.spectra_down[i, j, :]
        spectrum_atom_x = self.spectra_x[i, j, :]
            
        '''Asign resulting spectra'''      
        self.spectrum_atom = spectrum_atom        
        self.spectrum_atom_up = spectrum_atom_up     
        self.spectrum_atom_down = spectrum_atom_down
        self.spectrum_atom_x = spectrum_atom_x
            
        #select a fixed energy
        ####find peaks
        ndexes = dp.detect_peaks(spectrum_atom)#find the peaks
        peaks = self.vv[ndexes]  
        self.peaks = peaks
    
        if (len(peaks) == 0):
            i = 0
            ndexes = [0]
            energy = 10
        else:
            minpeak = min(abs(peaks))#find the minimum
            peaks2 = abs(peaks)
            peaks = peaks.tolist()
            peaks2 = peaks2.tolist()
            #i=peaks2.index(abs(minpeak))
            ii = find_nearest(peaks, -abs(minpeak))
            #print(i)
            #i=peaks2.index(minpeak) + 1#the index of the peak closest to zero
            energy = self.vv[ndexes[ii]]
            
        self.energy = energy
        self.ind = ndexes[ii]
        
        profile = self.spectra[i,:,ndexes[ii]]#profile at energy          
        self.profile = profile        
        
        
    def maps_2D(self, indx):        
        z = np.zeros([self.N_y, self.N_x], dtype = float)
        z_z = np.zeros([self.N_y, self.N_x], dtype = float)
        z_x = np.zeros([self.N_y, self.N_x], dtype = float)
        
        spectra_i = self.spectra[:,:, indx]
        spectra_i_up = self.spectra_up[:,:, indx]
        spectra_i_dn = self.spectra_down[:,:, indx]
        spectra_i_x = self.spectra_x[:,:, indx]        
        
        for j_atom in range(self.N_y):
            for i_atom in range(self.N_x):         
                
                z[j_atom,i_atom] = spectra_i[j_atom,i_atom]
                z_z[j_atom,i_atom] = spectra_i_up[j_atom,i_atom] - spectra_i_dn[j_atom,i_atom]
                z_x[j_atom,i_atom] = spectra_i_x[j_atom,i_atom]  
                
        self.z = z
        self.z_z = z_z
        self.z_x = z_x
                
        e =  np.zeros([self.N_x, self.N_omega], dtype = float)
        e_z =  np.zeros([self.N_x, self.N_omega], dtype = float)        
        e_x =  np.zeros([self.N_x, self.N_omega], dtype = float)
        
        ####selct spectra along x    
        spectra_i_x = self.spectra[self.row,:,:]
        spectra_i_xu = self.spectra_up[self.row,:,:]
        spectra_i_xd = self.spectra_down[self.row,:,:]
        spectra_i_xx = self.spectra_x[self.row,:,:]
        
        for i_omega in range(self.N_omega):
            for i_atom in range(self.N_x):
                e[i_atom,i_omega] = spectra_i_x[i_atom,i_omega]    
                e_z[i_atom,i_omega] = spectra_i_xu[i_atom,i_omega] - spectra_i_xd[i_atom,i_omega]    
                e_x[i_atom,i_omega] = spectra_i_xx[i_atom,i_omega]
                
        self.e = e
        self.e_z = e_z
        self.e_x = e_x























