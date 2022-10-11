#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 13:06:06 2022

@author: cristina
"""

#import matplotlib as mpl
#mpl.use('Agg') ####to run in Oberon/other clusters

import Spin_chain_class as clas

'''init class with total G calculation'''
Chain = clas.Shiba_chain()


'''Calculate spectra everywere in the array'''
Chain.spectra_calc()

'''Select spectrum we want to plot'''#eg. [i,j] = [row, borde]
Chain.select_spectra(Chain.row, Chain.borde)

#print("{:.8f}".format(Chain.peaks))
pp = Chain.peaks

'''plot spectra'''
import plot_spectra as ps
ps.plot_spectra(Chain.vv, Chain.spectrum_atom, Chain.spectrum_atom_up, Chain.spectrum_atom_down,\
                Chain.profile, Chain.energy)

'''create 2D plots'''
Chain.maps_2D(Chain.ind)

'''2-D plots'''
import plot_2D as p2
p2.plot_2D(Chain.z, Chain.z_z, Chain.z_x, Chain.energy, Chain.e, Chain.e_z, Chain.e_x, Chain.vv, Chain.N_x)

#print(Chain.peaks)
