#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#array dimensions
d = 1.0#distance between sites MUST BE 1.0!!!!!!
N_atoms = 25#number of atoms
borde = 3
wide = 3
layers = 1
lattice_param = 3.36

#parameters
alpha = 3.0#SOC eVA
state = 'FM'#spin state FM, AFM or helix
N_period = 3#spin helix 2=AFM
k_F = 0.183#Fermi wave vector
U = -2000./27211.6#%potential scatt#
j = 1700./27211.6#magnetic coupling
DOS = 1.0#normal density of states
s = 5.0/2.0#spin
delta = 0.75/27211.6 #SC gap
N_omega = 1001#numer of energy points
omega0 = -1.0#intial
omega1 = 1.0#final
Dynes = 0.01#Dynes parameter